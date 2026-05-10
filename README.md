# GLP-1R Binder Design — OMTX Hub Hackathon 2026

> **Track A · Biologics** — Automated peptide/protein binder design pipeline against GLP-1R, a validated drug target for type 2 diabetes and obesity.

Built for the [OMTX Hub](https://omtx.ai) GLP-1R Drug Design Hackathon. The pipeline generates semaglutide-derived peptide variants, submits them to cloud ML models, and tracks results against the live leaderboard.

---

## Results

### Leaderboard (Hub → BoltzGen tab)

| Rank | Variant | Model | iPTM | Notes |
|------|---------|-------|------|-------|
| 🥇 | `short_25_35_v2` | BoltzGen | **0.794** | Short peptide (25–35 aa), second sample |
| 🥈 | `short_20_30_v2` | BoltzGen | 0.731 | GLP-1 analog length range |
| 🥉 | `long_65_85` | BoltzGen | 0.549 | Miniprotein range |
| 4 | `long_75_90` | BoltzGen | 0.543 | Full miniprotein |
| 5 | `long_55_75` | BoltzGen | 0.431 | Long binder |

### AlphaFold2 Screening (peptide variant ranking, not on leaderboard)

| Variant | iPTM | Mutation | Insight |
|---------|------|----------|---------|
| T7A | **0.84** | Thr→Ala at pos 7 | Top hit — removes polar group at N-terminus |
| R29K | **0.84** | Arg→Lys at pos 29 | Conservative charge-preserving swap |
| E22D | 0.83 | Glu→Asp at pos 22 | Shorter side chain improves fit |
| Y13L | 0.82 | Tyr→Leu at pos 13 | Hydrophobic swap at helix core |
| A18C | 0.81 | Ala→Cys at pos 18 | Potential disulfide anchor site |

> iPTM (interface predicted TM-score) measures predicted binding quality. Values above 0.6 are strong; 0.8+ is excellent.

---

## How It Works

```
Semaglutide baseline (31 aa)
        │
        ▼
 variants.py ──► 117 peptide variants
   (single-point mutations, hotspot combos,
    N-term truncations, C-term extensions)
        │
        ▼
  run.py CLI
   ├─ submit-batch  ──► AlphaFold2 multimer  ($3/run) ── iPTM score
   ├─ boltzgen-batch ──► BoltzGen            ($5/run) ── iPTM → Leaderboard ✓
   └─ bindcraft-batch ─► BindCraft           ($10/run) ─ confidence → Leaderboard ✓
        │
        ▼
 results_tracker.py ──► results/results.json
        │
        ▼
  run.py leaderboard  (local ranked table)
```

**Two tracks on the Hub leaderboard:**
- **BoltzGen** — generative binder design from receptor CIF, varied length ranges
- **BindCraft** — de novo binder design from receptor PDB + hotspot residues

AlphaFold2 is used as a cheap screening tool to rank semaglutide variants before committing budget to the leaderboard models.

---

## Setup

```bash
git clone https://github.com/abhay-lal/glp1r-binder-design
cd glp1r-binder-design

python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

Download structure files (see [`data/README.md`](data/README.md)):
```bash
mkdir -p Hackathon_Files/Protein_Binder_Hackathon_Files
curl -o Hackathon_Files/Protein_Binder_Hackathon_Files/7ki0.cif \
  "https://files.rcsb.org/download/7KI0.cif"
curl -o Hackathon_Files/Protein_Binder_Hackathon_Files/7KI0_GLP1R_chainR.pdb \
  "https://files.rcsb.org/download/7KI0.pdb"
```

Set your OMTX API key (get it from Hub → User Menu → API Keys):
```bash
export OMTX_API_KEY=your_key_here
```

---

## Usage

```bash
# Check wallet balance
python pipeline/run.py credits

# Preview all 117 generated variants
python pipeline/run.py variants

# ── LEADERBOARD COMMANDS ────────────────────────────────────────
# BoltzGen: generative binder design ($5/run) → Hub leaderboard
python pipeline/run.py boltzgen-batch --limit 5

# BindCraft: de novo binder from hotspot residues ($10/run) → Hub leaderboard
python pipeline/run.py bindcraft-batch --limit 3

# Dry run to preview without spending
python pipeline/run.py boltzgen-batch --dry-run
python pipeline/run.py bindcraft-batch --dry-run

# ── SCREENING COMMANDS ──────────────────────────────────────────
# AlphaFold2 multimer: cheap variant ranking ($3/run, not on leaderboard)
python pipeline/run.py submit-batch --strategy hotspot --model alphafold --limit 7

# Chai-1: even cheaper screening ($0.85/run)
python pipeline/run.py submit-batch --strategy hotspot --model chai1 --limit 7

# ── RESULTS ────────────────────────────────────────────────────
python pipeline/run.py poll          # update job statuses + backfill scores
python pipeline/run.py leaderboard   # ranked local results table
python pipeline/run.py best --top 5  # top 5 by score
python pipeline/run.py plan          # budget breakdown + recommendations
```

---

## Pipeline Modules

| File | Purpose |
|------|---------|
| `pipeline/sequence.py` | GLP-1R sequence (384 aa), semaglutide baseline, orthosteric pocket contacts |
| `pipeline/variants.py` | 117 semaglutide variant generator (mutations, truncations, extensions) |
| `pipeline/leaderboard_jobs.py` | BindCraft hotspot combinations + BoltzGen length sweep configs |
| `pipeline/api_client.py` | OMTX Hub SDK wrapper — submit, poll, score extraction for all models |
| `pipeline/results_tracker.py` | Local results store with deduplication and leaderboard ranking |
| `pipeline/yaml_builder.py` | Boltz-2 YAML reference files (for manual inspection) |
| `pipeline/run.py` | CLI entrypoint — all commands |

---

## Variant Design Space

| Strategy | Count | Description |
|----------|-------|-------------|
| `single_point` | 74 | Conservative residue swaps across all 31 positions |
| `double_hotspot` | 19 | Combinations of top single-point mutations |
| `nterm_truncation` | 11 | N-terminal trimming (test minimal binding fragment) |
| `hotspot` | 7 | Literature-guided mutations at GLP-1R contact residues |
| `cterm_extension` | 5 | Helicity-enhancing C-terminal extensions |
| `baseline` | 1 | Unmodified semaglutide analog |
| **Total** | **117** | |

---

## Target Information

- **Receptor:** GLP-1R (PDB [7KI0](https://www.rcsb.org/structure/7KI0)) — Class B GPCR, 384 resolved residues
- **Baseline peptide:** `HAEGTFTSDVSSYLEGQAAKEFIAWLVRGRG` (31 aa semaglutide analog from 7KI0)
- **Key binding pocket residues:** Y152, R190, Y241, W306, R310, L314, H363, E364, L379, R380, F381, L388
- **BoltzGen chain:** `F` in `7ki0.cif` (entity 6, label_asym_id F — full RCSB mmCIF required)
- **BindCraft chain:** `R` in `7KI0_GLP1R_chainR.pdb` (pre-extracted receptor)

---

## Model Pricing (OMTX Hub, 2026)

| Model | Cost | Leaderboard | Use case |
|-------|------|-------------|----------|
| BoltzGen | $5.00/run | ✅ Track A | Generative binder design |
| BindCraft | $10.00/run | ✅ Track A | De novo binder from hotspots |
| AlphaFold2 | $3.00/run | ❌ screening | Multimer iPTM scoring |
| Chai-1 | $0.85/run | ❌ screening | Cheap multimer screening |
| Boltz-2 | $0.85/run | ✅ Track B | Protein-ligand ΔG |
