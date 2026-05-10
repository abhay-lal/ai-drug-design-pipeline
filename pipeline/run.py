#!/usr/bin/env python3
"""
GLP-1R Track A Pipeline — OMTX Hub Hackathon 2026

LEADERBOARD MODELS (these post scores to Hub → Leaderboard):
  bindcraft  — De novo binder design from receptor PDB + hotspot residues
               Vary hotspot subsets → each job designs a novel binder
               Score: confidence — HIGHER is better  ← Track A leaderboard
  boltzgen   — Generative binder design from receptor CIF
               Vary binder length range → each job generates new binders
               Score: confidence — HIGHER is better  ← Track A leaderboard

SCREENING MODELS (useful for scoring your peptide variants, not on leaderboard):
  alphafold  — AlphaFold multimer: receptor + your peptide → pTM/ipTM score
               Use to rank your semaglutide mutations before guiding BindCraft

Usage:
  # Check wallet balance
  python run.py credits

  # ── LEADERBOARD COMMANDS ──────────────────────────────────────────────
  # Submit all BindCraft hotspot combinations (hits leaderboard)
  python run.py bindcraft-batch --limit 10 --wait

  # Submit all BoltzGen length sweeps (hits leaderboard)
  python run.py boltzgen-batch --limit 10 --wait

  # Preview what jobs would be submitted
  python run.py bindcraft-batch --dry-run
  python run.py boltzgen-batch --dry-run

  # ── SCREENING COMMANDS ────────────────────────────────────────────────
  # Score your peptide variants with AlphaFold multimer
  python run.py submit-batch --strategy hotspot --model alphafold --limit 7 --wait
  python run.py variants   # list all 117 peptide variants

  # ── RESULTS ───────────────────────────────────────────────────────────
  python run.py poll
  python run.py leaderboard
  python run.py best --top 5

Set your API key:
  export OMTX_API_KEY=your_key_here
  Get it: Hub → User Menu → API Keys → omtx.ai
"""

import os
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import click
from rich.console import Console
from rich.table import Table

sys.path.insert(0, str(Path(__file__).parent))

from api_client import OmtxClient, PRICES, LEADERBOARD_MODELS
from sequence import GLP1R_SEQUENCE
from variants import generate_all_variants, PeptideVariant, BASELINE
from yaml_builder import write_yaml_job, write_all_yaml_jobs
from results_tracker import record_result, load_leaderboard, print_leaderboard, get_best, already_submitted, submitted_names
from leaderboard_jobs import generate_bindcraft_jobs, generate_boltzgen_jobs, BindCraftJob, BoltzGenJob


JOBS_DIR    = Path(__file__).parent.parent / "jobs"
RESULTS_DIR = Path(__file__).parent.parent / "results"
PDB_PATH    = (
    Path(__file__).parent.parent
    / "Hackathon_Files/Protein_Binder_Hackathon_Files/7KI0_GLP1R_chainR.pdb"
)
CIF_PATH    = (
    Path(__file__).parent.parent
    / "Hackathon_Files/Protein_Binder_Hackathon_Files/7ki0.cif"
)
# GLP-1R is chain F in the original 7ki0.cif (entity 6, label_asym_id F)
# The chainR.cif files are stripped exports that lack _entity_poly_seq required by BoltzGen
BOLTZGEN_CHAIN = "F"

console = Console()

# Cache uploaded artifact IDs so we don't re-upload on every batch call
_artifact_cache: dict[str, str] = {}


# ── Helpers ──────────────────────────────────────────────────────────────────

def _get_client():
    from api_client import OmtxClient
    key = os.environ.get("OMTX_API_KEY", "")
    if not key:
        console.print("[red bold]OMTX_API_KEY not set.[/red bold]")
        console.print("  [dim]export OMTX_API_KEY=your_key_here[/dim]")
        console.print("  [dim]Get it: Hub → User Menu → API Keys → omtx.ai[/dim]")
        sys.exit(1)
    return OmtxClient(api_key=key)


def _ensure_pdb_artifact(client) -> str:
    """Upload receptor PDB once and cache the artifact_id."""
    if "pdb" not in _artifact_cache:
        console.print(f"Uploading receptor PDB: [cyan]{PDB_PATH.name}[/cyan] …")
        _artifact_cache["pdb"] = client.upload_pdb(PDB_PATH)
        console.print(f"  artifact_id: [green]{_artifact_cache['pdb']}[/green]")
    return _artifact_cache["pdb"]


def _ensure_cif_artifact(client) -> str:
    """Upload the full 7ki0.cif (from RCSB) once and cache the artifact_id.

    BoltzGen requires a complete mmCIF with _entity_poly_seq.entity_id and
    other standard tables.  The stripped chainR.cif exports lack these.
    GLP-1R is chain F (entity 6) in the original 7ki0.cif.
    """
    if "cif" not in _artifact_cache:
        console.print(f"Uploading full receptor CIF: [cyan]{CIF_PATH.name}[/cyan] …")
        _artifact_cache["cif"] = client.upload_cif(CIF_PATH)
        console.print(f"  artifact_id: [green]{_artifact_cache['cif']}[/green]")
    return _artifact_cache["cif"]


def _variants_table(variants: list[PeptideVariant]) -> Table:
    t = Table(title=f"Peptide Variants ({len(variants)} total)", show_lines=False)
    t.add_column("#",        style="dim",    width=5)
    t.add_column("Name",     style="cyan",   width=22)
    t.add_column("Strategy", style="yellow", width=18)
    t.add_column("Len",      justify="right",width=5)
    t.add_column("Sequence", style="green",  width=50)
    for i, v in enumerate(variants):
        t.add_row(str(i), v.name, v.strategy, str(v.length), v.sequence)
    return t


def _submit_one(client, variant: PeptideVariant, model: str, pdb_artifact_id: str | None) -> "JobResult":
    if model == "chai1":
        # Cheap complex scoring: $0.85/run vs AlphaFold $3.00/run
        return client.submit_chai1(
            receptor_sequence=GLP1R_SEQUENCE,
            binder_sequence=variant.sequence,
            variant_name=variant.name,
        )
    elif model == "alphafold":
        # AlphaFold multimer: $3.00/run (use chai1 instead to save cost)
        return client.submit_alphafold(
            receptor_sequence=GLP1R_SEQUENCE,
            binder_sequence=variant.sequence,
            variant_name=variant.name,
        )
    elif model == "bindcraft":
        return client.submit_bindcraft(
            receptor_pdb_path=PDB_PATH,
            variant_name=variant.name,
            artifact_id=pdb_artifact_id,
        )
    elif model == "boltzgen":
        cif_id = _ensure_cif_artifact(client)
        return client.submit_boltzgen(
            target_cif_artifact_id=cif_id,
            variant_name=variant.name,
            target_chain_id=BOLTZGEN_CHAIN,
        )
    elif model == "boltz2":
        return client.submit_boltz2(
            protein_sequence=GLP1R_SEQUENCE,
            variant_name=variant.name,
        )
    else:
        raise ValueError(f"Unknown model: {model}. Choose: chai1 | alphafold | bindcraft | boltzgen | boltz2")


# ── Commands ──────────────────────────────────────────────────────────────────

@click.group()
def cli():
    """GLP-1R Track A — OMTX Hub Hackathon 2026."""


@cli.command()
@click.option("--conservative/--all-subs", default=True)
@click.option("--no-combos", is_flag=True)
def variants(conservative, no_combos):
    """Preview all generated variants."""
    all_v = generate_all_variants(
        conservative_only=conservative,
        include_combos=not no_combos,
    )
    console.print(_variants_table(all_v))
    console.print(f"\n[bold]Baseline:[/bold] {BASELINE} ({len(BASELINE)} aa)")
    console.print(f"[bold]Total:[/bold] {len(all_v)} variants\n")
    strats = {}
    for v in all_v:
        strats[v.strategy] = strats.get(v.strategy, 0) + 1
    for s, n in sorted(strats.items()):
        console.print(f"  {s:<22} {n}")


@cli.command()
@click.option("--conservative/--all-subs", default=True)
@click.option("--no-combos", is_flag=True)
@click.option("--output-dir", default=str(JOBS_DIR), show_default=True)
def generate(conservative, no_combos, output_dir):
    """Write YAML reference files for all variants (for inspection/manual upload)."""
    all_v = generate_all_variants(
        conservative_only=conservative,
        include_combos=not no_combos,
    )
    out = Path(output_dir)
    console.print(f"Writing {len(all_v)} YAML files to [cyan]{out}[/cyan] …")
    paths = write_all_yaml_jobs(all_v, out)
    console.print(f"[green]Done.[/green] {len(paths)} files written.")


@cli.command()
def credits():
    """Check your OMTX wallet balance."""
    client = _get_client()
    cents = client.check_credits()
    dollars = cents / 100
    console.print(f"[bold]Wallet balance:[/bold] [green]${dollars:.2f}[/green] ({int(cents)} credits)")


@cli.command()
@click.option("--budget", default=50.0, show_default=True, help="Total budget in USD")
def plan(budget):
    """Show pricing, budget status, and optimal job allocation."""
    import json

    # Tally spend from results.json
    results_file = RESULTS_DIR / "results.json"
    spend_by_model: dict[str, float] = {}
    count_by_model: dict[str, int] = {}
    if results_file.exists():
        with open(results_file) as f:
            for r in json.load(f):
                m = r.get("model", "unknown")
                p = PRICES.get(m, 0)
                spend_by_model[m] = spend_by_model.get(m, 0) + p
                count_by_model[m] = count_by_model.get(m, 0) + 1

    total_spent = sum(spend_by_model.values())
    remaining   = budget - total_spent

    # Summary stats
    console.print(f"\n[bold]Budget:[/bold]    ${budget:.2f}")
    console.print(f"[bold]Spent:[/bold]     [red]${total_spent:.2f}[/red]  ({sum(count_by_model.values())} jobs)")
    console.print(f"[bold]Remaining:[/bold] [green]${remaining:.2f}[/green]\n")

    # Spend breakdown
    if spend_by_model:
        t = Table(title="Spend by Model", show_lines=False)
        t.add_column("Model",  style="cyan",   width=16)
        t.add_column("Jobs",   justify="right",width=6)
        t.add_column("$/run",  justify="right",width=8)
        t.add_column("Total",  justify="right",width=10, style="yellow")
        t.add_column("Leaderboard", width=30, style="dim")
        for m, total in sorted(spend_by_model.items(), key=lambda x: -x[1]):
            lb = LEADERBOARD_MODELS.get(m, "—  (screening only)")
            t.add_row(m, str(count_by_model[m]), f"${PRICES.get(m,0):.2f}", f"${total:.2f}", lb)
        console.print(t)

    # What fits in remaining budget
    console.print(f"\n[bold]With ${remaining:.2f} remaining you can run:[/bold]")
    t2 = Table(show_lines=False, show_header=False)
    t2.add_column("Model",  style="cyan",  width=16)
    t2.add_column("Runs",   justify="right", width=6)
    t2.add_column("Cost",   justify="right", width=10, style="yellow")
    t2.add_column("Note",   width=40, style="dim")
    priority = [
        ("bindcraft",   "LEADERBOARD Track A — de novo design"),
        ("boltzgen",    "LEADERBOARD Track A — generative design"),
        ("chai1",       "Cheap peptide screening (vs alphafold $3)"),
        ("openfold3",   "Multimer complex scoring"),
        ("boltz2",      "LEADERBOARD Track B — ligand ΔG"),
        ("alphafold",   "Multimer (use chai1 instead to save cost)"),
    ]
    for m, note in priority:
        price = PRICES.get(m, 0)
        n = int(remaining // price) if price else 0
        t2.add_row(m, str(n), f"${price:.2f}/run", note)
    console.print(t2)

    # Recommended allocation
    console.print(f"\n[bold]Recommended allocation (${remaining:.2f}):[/bold]")
    bc = min(2, int(remaining // PRICES["bindcraft"]))
    left = remaining - bc * PRICES["bindcraft"]
    bg = min(3, int(left // PRICES["boltzgen"]))
    left -= bg * PRICES["boltzgen"]
    console.print(f"  [cyan]bindcraft-batch --limit {bc}[/cyan]   ${bc * PRICES['bindcraft']:.2f}  → leaderboard")
    console.print(f"  [cyan]boltzgen-batch  --limit {bg}[/cyan]   ${bg * PRICES['boltzgen']:.2f}  → leaderboard")
    console.print(f"  [dim]Buffer: ${left:.2f}[/dim]")


@cli.command()
@click.option("--variant", default="baseline", show_default=True,
              help="Variant name — run 'variants' to list all names")
@click.option("--model", default="chai1", show_default=True,
              help="Model: chai1 ($0.85) | alphafold ($3.00) | bindcraft ($10) | boltzgen ($5) | boltz2 ($0.85)")
@click.option("--wait/--no-wait", default=False)
def submit(variant, model, wait):
    """Submit a single variant."""
    client = _get_client()
    all_v = generate_all_variants()
    match = next((v for v in all_v if v.name == variant), None)
    if not match:
        console.print(f"[red]Unknown variant:[/red] {variant}")
        console.print("Run [cyan]python run.py variants[/cyan] to see all names.")
        sys.exit(1)

    pdb_artifact_id = _ensure_pdb_artifact(client) if model == "bindcraft" else None

    console.print(f"Submitting [cyan]{match.name}[/cyan] ({match.length} aa) → [yellow]{model}[/yellow] …")
    result = _submit_one(client, match, model, pdb_artifact_id)
    record_result(result)
    console.print(f"[green]Queued![/green] Job ID: [bold]{result.job_id}[/bold]")

    if wait:
        console.print("Waiting for result …")
        result = client.poll_until_done(result)
        record_result(result)
        status_color = "green" if result.status.value == "completed" else "red"
        console.print(f"[{status_color}]{result}[/{status_color}]")


@cli.command("submit-batch")
@click.option("--strategy", default=None,
              help="Filter: single_point | cterm_extension | nterm_truncation | hotspot | double_hotspot | baseline")
@click.option("--model", default="chai1", show_default=True)
@click.option("--limit", default=20, show_default=True)
@click.option("--workers", default=4, show_default=True)
@click.option("--wait/--no-wait", default=False)
def submit_batch(strategy, model, limit, workers, wait):
    """Submit multiple variants concurrently."""
    client = _get_client()
    all_v = generate_all_variants()
    if strategy:
        all_v = [v for v in all_v if v.strategy == strategy]
    all_v = all_v[:limit]

    # Pre-upload PDB once if needed
    pdb_artifact_id = _ensure_pdb_artifact(client) if model == "bindcraft" else None

    console.print(
        f"Submitting [bold]{len(all_v)}[/bold] variants "
        f"(strategy=[cyan]{strategy or 'all'}[/cyan], model=[yellow]{model}[/yellow]) …"
    )

    submitted = []

    def _worker(v: PeptideVariant):
        result = _submit_one(client, v, model, pdb_artifact_id)
        record_result(result)
        return result

    with ThreadPoolExecutor(max_workers=workers) as ex:
        futures = {ex.submit(_worker, v): v for v in all_v}
        for future in as_completed(futures):
            v = futures[future]
            try:
                result = future.result()
                submitted.append(result)
                console.print(f"  [green]✓[/green] {v.name:<22} → {result.job_id}")
            except Exception as e:
                console.print(f"  [red]✗[/red] {v.name}: {e}")

    console.print(f"\n[bold]{len(submitted)}/{len(all_v)}[/bold] jobs queued.")

    if wait and submitted:
        console.print("\nPolling until all done …")
        for result in submitted:
            result = client.poll_until_done(result, poll_interval=20)
            record_result(result)
            color = "green" if result.status.value == "completed" else "red"
            console.print(f"  [{color}]{result}[/{color}]")


@cli.command()
@click.option("--interval", default=20, show_default=True)
def poll(interval):
    """Poll all pending/running jobs and update local results."""
    import json
    results_file = RESULTS_DIR / "results.json"
    if not results_file.exists():
        console.print("No results file. Submit some jobs first.")
        return

    with open(results_file) as f:
        all_results = json.load(f)

    pending = [
        r for r in all_results
        if r.get("status") in ("pending", "queued", "running", "unknown")
    ]
    # Also backfill completed jobs that are missing a score
    scoreless = [
        r for r in all_results
        if r.get("status") == "completed" and r.get("score") is None
    ]

    if not pending and not scoreless:
        console.print("No pending jobs and all scores are populated.")
        return

    client = _get_client()

    if pending:
        console.print(f"Polling [bold]{len(pending)}[/bold] in-flight jobs …")
        for r in pending:
            try:
                updated = client.get_job(
                    r["job_id"],
                    r.get("variant_name", ""),
                    r.get("model", ""),
                )
                record_result(updated)
                color = "green" if updated.status.value == "completed" else "yellow"
                score_str = f"  score={updated.score:.4f}" if updated.score else ""
                console.print(
                    f"  [{color}]{updated.status.value:<10}[/{color}] "
                    f"{updated.variant_name:<22}{score_str}"
                )
            except Exception as e:
                console.print(f"  [red]Error polling {r['job_id']}:[/red] {e}")
            time.sleep(0.3)

    if scoreless:
        console.print(f"\nBackfilling scores for [bold]{len(scoreless)}[/bold] completed jobs …")
        for r in scoreless:
            try:
                updated = client.get_job(
                    r["job_id"],
                    r.get("variant_name", ""),
                    r.get("model", ""),
                )
                record_result(updated)
                score_str = f"  score={updated.score:.4f} [{updated.score_label}]" if updated.score else "  no score in payload"
                console.print(f"  [cyan]{r.get('variant_name','?'):<22}[/cyan] ({r.get('model')}){score_str}")
            except Exception as e:
                console.print(f"  [red]Error backfilling {r['job_id']}:[/red] {e}")
            time.sleep(0.3)


@cli.command()
@click.option("--model", default=None)
def leaderboard(model):
    """Print local leaderboard of completed jobs."""
    print_leaderboard(model)


@cli.command()
@click.option("--top", default=5, show_default=True)
@click.option("--model", default=None)
def best(top, model):
    """Show top-N best-scoring variants."""
    board = load_leaderboard(model)[:top]
    if not board:
        console.print("No completed results yet.")
        return

    t = Table(title=f"Top {top} Variants", show_lines=False)
    t.add_column("Rank",    width=5)
    t.add_column("Variant", style="cyan",   width=25)
    t.add_column("Model",   style="yellow", width=12)
    t.add_column("Score",   justify="right",width=10)
    t.add_column("Metric",  width=12)
    t.add_column("Job ID",  style="dim",    width=36)

    for rank, r in enumerate(board, 1):
        score_str = f"{r['score']:.4f}" if r.get("score") is not None else "—"
        t.add_row(
            str(rank),
            r.get("variant_name", "?"),
            r.get("model", "?"),
            score_str,
            r.get("score_label", "?"),
            r.get("job_id", "?"),
        )
    console.print(t)


# ── Leaderboard commands ──────────────────────────────────────────────────────

@cli.command("bindcraft-batch")
@click.option("--limit", default=1, show_default=True,
              help="Max new BindCraft jobs to submit ($10 each)")
@click.option("--workers", default=2, show_default=True)
@click.option("--wait/--no-wait", default=False)
@click.option("--dry-run", is_flag=True, help="Preview jobs without submitting")
def bindcraft_batch(limit, workers, wait, dry_run):
    """
    Submit BindCraft de novo design jobs — posts scores to Hub → Leaderboard.
    Cost: $10/run. Skips any job name already in results.json.
    """
    from api_client import OmtxClient
    all_jobs = generate_bindcraft_jobs()

    # Dedup: skip names already submitted for bindcraft
    done = submitted_names("bindcraft")
    new_jobs = [j for j in all_jobs if j.name not in done]
    skipped = len(all_jobs) - len(new_jobs)
    new_jobs = new_jobs[:limit]

    cost = len(new_jobs) * 10
    t = Table(title=f"BindCraft Jobs — {len(new_jobs)} new  (${cost})  [{skipped} already submitted, skipped]", show_lines=False)
    t.add_column("#",        width=4,  style="dim")
    t.add_column("Name",     width=22, style="cyan")
    t.add_column("Hotspots", width=50, style="green")
    t.add_column("Cost",     width=6,  justify="right", style="yellow")
    for i, j in enumerate(new_jobs):
        t.add_row(str(i+1), j.name, j.hotspot_residues, "$10")
    console.print(t)
    console.print(f"[bold]Estimated cost: [yellow]${cost}[/yellow][/bold]")

    if dry_run:
        console.print("[dim]Dry run — no jobs submitted.[/dim]")
        return

    if not new_jobs:
        console.print("[dim]Nothing new to submit.[/dim]")
        return

    client = _get_client()
    pdb_id = _ensure_pdb_artifact(client)

    submitted = []

    def _submit(j: BindCraftJob):
        result = client.submit_bindcraft(
            receptor_pdb_path=PDB_PATH,
            variant_name=j.name,
            hotspot_residues=j.hotspot_residues,
            target_chains="R",
            artifact_id=pdb_id,
        )
        record_result(result)
        return result

    console.print(f"\nSubmitting [bold]{len(new_jobs)}[/bold] BindCraft jobs …")
    with ThreadPoolExecutor(max_workers=workers) as ex:
        futures = {ex.submit(_submit, j): j for j in new_jobs}
        for future in as_completed(futures):
            j = futures[future]
            try:
                result = future.result()
                submitted.append(result)
                console.print(f"  [green]✓[/green] {j.name:<22} → {result.job_id}")
            except Exception as e:
                console.print(f"  [red]✗[/red] {j.name}: {e}")

    console.print(f"\n[bold]{len(submitted)}/{len(new_jobs)}[/bold] BindCraft jobs queued.")
    console.print("[dim]Scores → Hub → Leaderboard → BindCraft tab.[/dim]")

    if wait and submitted:
        console.print("Polling until done (BindCraft takes ~10-20 min) …")
        for result in submitted:
            result = client.poll_until_done(result, poll_interval=30, timeout=2400)
            record_result(result)
            color = "green" if result.status.value == "completed" else "red"
            console.print(f"  [{color}]{result}[/{color}]")


@cli.command("boltzgen-batch")
@click.option("--limit", default=3, show_default=True,
              help="Max new BoltzGen jobs to submit ($5 each)")
@click.option("--workers", default=3, show_default=True)
@click.option("--wait/--no-wait", default=False)
@click.option("--dry-run", is_flag=True)
def boltzgen_batch(limit, workers, wait, dry_run):
    """
    Submit BoltzGen generative design jobs — posts scores to Hub → Leaderboard.
    Cost: $5/run. Skips any job name already in results.json.
    """
    all_jobs = generate_boltzgen_jobs()

    done = submitted_names("boltzgen")
    new_jobs = [j for j in all_jobs if j.name not in done]
    skipped = len(all_jobs) - len(new_jobs)
    new_jobs = new_jobs[:limit]

    cost = len(new_jobs) * 5
    t = Table(title=f"BoltzGen Jobs — {len(new_jobs)} new  (${cost})  [{skipped} already submitted, skipped]", show_lines=False)
    t.add_column("#",        width=4,  style="dim")
    t.add_column("Name",     width=22, style="cyan")
    t.add_column("Min len",  width=8,  justify="right")
    t.add_column("Max len",  width=8,  justify="right")
    t.add_column("Cost",     width=6,  justify="right", style="yellow")
    for i, j in enumerate(new_jobs):
        t.add_row(str(i+1), j.name, str(j.binder_length_min), str(j.binder_length_max), "$5")
    console.print(t)
    console.print(f"[bold]Estimated cost: [yellow]${cost}[/yellow][/bold]")

    if dry_run:
        console.print("[dim]Dry run — no jobs submitted.[/dim]")
        return

    if not new_jobs:
        console.print("[dim]Nothing new to submit.[/dim]")
        return

    client = _get_client()
    cif_id = _ensure_cif_artifact(client)

    submitted = []

    def _submit(j: BoltzGenJob):
        result = client.submit_boltzgen(
            target_cif_artifact_id=cif_id,
            variant_name=j.name,
            binder_length_min=j.binder_length_min,
            binder_length_max=j.binder_length_max,
            target_chain_id=BOLTZGEN_CHAIN,
        )
        record_result(result)
        return result

    console.print(f"\nSubmitting [bold]{len(new_jobs)}[/bold] BoltzGen jobs …")
    with ThreadPoolExecutor(max_workers=workers) as ex:
        futures = {ex.submit(_submit, j): j for j in new_jobs}
        for future in as_completed(futures):
            j = futures[future]
            try:
                result = future.result()
                submitted.append(result)
                console.print(f"  [green]✓[/green] {j.name:<22} → {result.job_id}")
            except Exception as e:
                console.print(f"  [red]✗[/red] {j.name}: {e}")

    console.print(f"\n[bold]{len(submitted)}/{len(new_jobs)}[/bold] BoltzGen jobs queued.")
    console.print("[dim]Scores → Hub → Leaderboard → BoltzGen tab.[/dim]")

    if wait and submitted:
        console.print("Polling until done …")
        for result in submitted:
            result = client.poll_until_done(result, poll_interval=30, timeout=2400)
            record_result(result)
            color = "green" if result.status.value == "completed" else "red"
            console.print(f"  [{color}]{result}[/{color}]")


if __name__ == "__main__":
    cli()
