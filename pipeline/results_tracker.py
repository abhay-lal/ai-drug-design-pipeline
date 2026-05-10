"""
Track submission results and maintain a local leaderboard.

Results are persisted to results/results.json so you can resume
across sessions and compare variants.
"""

import json
from datetime import datetime
from pathlib import Path

from api_client import JobResult, JobStatus

RESULTS_FILE = Path(__file__).parent.parent / "results" / "results.json"


def _load_results() -> list[dict]:
    if RESULTS_FILE.exists():
        with open(RESULTS_FILE) as f:
            return json.load(f)
    return []


def _save_results(results: list[dict]) -> None:
    RESULTS_FILE.parent.mkdir(parents=True, exist_ok=True)
    with open(RESULTS_FILE, "w") as f:
        json.dump(results, f, indent=2, default=str)


def record_result(job: JobResult) -> None:
    """Append or update a job result in the local results store."""
    results = _load_results()
    # Update existing entry if job_id matches
    for entry in results:
        if entry.get("job_id") == job.job_id:
            entry.update(_job_to_dict(job))
            _save_results(results)
            return
    # New entry
    results.append(_job_to_dict(job))
    _save_results(results)


def load_leaderboard(model: str | None = None) -> list[dict]:
    """
    Return completed jobs sorted by score (desc for confidence, asc for dG).
    Optionally filter by model.
    """
    results = _load_results()
    completed = [r for r in results if r.get("status") == JobStatus.COMPLETED.value]
    if model:
        completed = [r for r in completed if r.get("model") == model]

    def sort_key(r):
        score = r.get("score")
        label = r.get("score_label", "")
        if score is None:
            return (1, 0)
        # Lower is better for dG
        if label in ("dg", "delta_g"):
            return (0, score)
        # Higher is better for everything else
        return (0, -score)

    return sorted(completed, key=sort_key)


def print_leaderboard(model: str | None = None) -> None:
    """Print a formatted leaderboard table to stdout."""
    board = load_leaderboard(model)
    if not board:
        print("No completed results yet.")
        return

    header = f"{'Rank':<5} {'Variant':<25} {'Model':<14} {'Score':>10} {'Label':<12} {'Job ID'}"
    print(header)
    print("-" * len(header))
    for rank, r in enumerate(board, 1):
        score_str = f"{r['score']:.4f}" if r.get("score") is not None else "—"
        print(
            f"{rank:<5} {r.get('variant_name', '?'):<25} {r.get('model', '?'):<14} "
            f"{score_str:>10} {(r.get('score_label') or '?'):<12} {r.get('job_id', '?')}"
        )


def get_best(model: str | None = None) -> dict | None:
    board = load_leaderboard(model)
    return board[0] if board else None


def already_submitted(variant_name: str, model: str) -> bool:
    """Return True if a job for this variant+model is already recorded (any status)."""
    for r in _load_results():
        if r.get("variant_name") == variant_name and r.get("model") == model:
            return True
    return False


def submitted_names(model: str) -> set[str]:
    """Return the set of variant names already submitted for a given model."""
    return {
        r["variant_name"]
        for r in _load_results()
        if r.get("model") == model and r.get("variant_name")
    }


def _job_to_dict(job: JobResult) -> dict:
    d = {
        "job_id": job.job_id,
        "variant_name": job.variant_name,
        "model": job.model,
        "status": job.status.value,
        "score": job.score,
        "score_label": job.score_label,
        "error": job.error,
        "recorded_at": datetime.utcnow().isoformat(),
    }
    # Persist BindCraft diagnostics if present
    raw_resp = (job.raw.get("response_payload") or {})
    raw_metrics = raw_resp.get("metrics") or {}
    if raw_metrics.get("outcome"):
        d["outcome"] = raw_metrics["outcome"]
    if raw_metrics.get("termination_reason"):
        d["termination_reason"] = raw_metrics["termination_reason"]
    if raw_metrics.get("suggested_next_action"):
        d["suggested_next_action"] = raw_metrics["suggested_next_action"]
    return d
