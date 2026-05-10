"""
Generate BindCraft and BoltzGen job configurations for Track A leaderboard submissions.

BindCraft: Takes receptor PDB + hotspot residues → designs novel binders de novo.
           We vary which subset of the 12 orthosteric pocket residues to target.

BoltzGen:  Takes receptor CIF → generates binders from scratch.
           We vary binder length range and protocol to explore diverse designs.

Both models post directly to Hub → Leaderboard (BindCraft tab / BoltzGen tab).
"""

import itertools
from dataclasses import dataclass, field

# ── Pocket residues ──────────────────────────────────────────────────────────
# Full orthosteric binding pocket from the hackathon guide.
# Residue numbers are as they appear in 7KI0_GLP1R_chainR.pdb (PDB numbering).
# BindCraft hotspot format: "R<chain><resnum>" e.g. "RA152"
# Chain in the uploaded PDB is R, but BindCraft expects it as the chain label
# present in the file. Chain R → we specify target_chains="R".

POCKET_RESIDUES = {
    "Y152": 152, "R190": 190, "Y241": 241, "R310": 310,
    "W306": 306, "L314": 314, "E364": 364, "H363": 363,
    "L379": 379, "R380": 380, "F381": 381, "L388": 388,
}

# Core contact residues from the YAML constraints (most important for binding)
CORE_RESIDUES = ["Y152", "W306", "E364", "R380"]

# Secondary shell residues
SHELL_RESIDUES = ["R190", "Y241", "R310", "L314", "H363", "L379", "F381", "L388"]


# ── BindCraft job configs ─────────────────────────────────────────────────────

@dataclass
class BindCraftJob:
    name: str
    hotspot_residues: str      # e.g. "R152,R306,R364,R380" (chain R prefix)
    target_chains: str = "R"
    description: str = ""
    tags: list[str] = field(default_factory=list)

    def __str__(self) -> str:
        return f"{self.name}: hotspots=[{self.hotspot_residues}]"


def _fmt_hotspots(residue_names: list[str], chain: str = "R") -> str:
    """Format residue names as BindCraft hotspot string: 'R152,R306,...'"""
    nums = [str(POCKET_RESIDUES[r]) for r in residue_names]
    return ",".join(f"{chain}{n}" for n in nums)


def generate_bindcraft_jobs() -> list[BindCraftJob]:
    """
    Generate diverse BindCraft job configurations by varying hotspot subsets.
    Strategy: core-only, core+shell pairs, full pocket, and literature-guided combos.
    """
    jobs: list[BindCraftJob] = []
    seen: set[str] = set()

    def _add(name: str, residues: list[str], desc: str, tags: list[str]) -> None:
        hs = _fmt_hotspots(residues)
        if hs in seen:
            return
        seen.add(hs)
        jobs.append(BindCraftJob(
            name=name,
            hotspot_residues=hs,
            description=desc,
            tags=tags,
        ))

    # 1. Core 4 residues (from YAML constraints)
    _add("core4", CORE_RESIDUES,
         "Core orthosteric contacts: Y152, W306, E364, R380",
         ["core"])

    # 2. All 12 pocket residues
    _add("full_pocket", list(POCKET_RESIDUES.keys()),
         "Full 12-residue orthosteric pocket",
         ["full"])

    # 3. Core + each shell residue pair
    for r in SHELL_RESIDUES:
        _add(f"core+{r}", CORE_RESIDUES + [r],
             f"Core 4 + {r}",
             ["core", "extended"])

    # 4. Core + all shell pairs (2 at a time)
    for r1, r2 in itertools.combinations(SHELL_RESIDUES, 2):
        name = f"core+{r1}+{r2}"
        _add(name, CORE_RESIDUES + [r1, r2],
             f"Core 4 + {r1} + {r2}",
             ["core", "combo"])
        if len(jobs) >= 30:   # cap to keep wallet usage reasonable
            break

    # 5. Top half of pocket (N-terminal side, TM helices 1-4)
    top_half = ["Y152", "R190", "Y241", "R310", "W306", "L314"]
    _add("top_half", top_half,
         "Top half of pocket (TM1-4 residues)",
         ["half"])

    # 6. Bottom half (C-terminal side, TM helices 5-7)
    bot_half = ["E364", "H363", "L379", "R380", "F381", "L388"]
    _add("bot_half", bot_half,
         "Bottom half of pocket (TM5-7 residues)",
         ["half"])

    # 7. Charged residues only (ionic contacts)
    charged = ["R190", "R310", "E364", "R380"]
    _add("charged", charged,
         "Charged residues only (ionic binding network)",
         ["charged"])

    # 8. Aromatic residues only (hydrophobic core)
    aromatic = ["Y152", "W306", "Y241", "F381"]
    _add("aromatic", aromatic,
         "Aromatic residues (hydrophobic binding core)",
         ["aromatic"])

    return jobs


# ── BoltzGen job configs ──────────────────────────────────────────────────────

@dataclass
class BoltzGenJob:
    name: str
    binder_length_min: int
    binder_length_max: int
    protocol: str = "protein_anything"
    target_chain_id: str = "R"
    description: str = ""
    tags: list[str] = field(default_factory=list)

    def __str__(self) -> str:
        return f"{self.name}: len=[{self.binder_length_min}-{self.binder_length_max}] protocol={self.protocol}"


def generate_boltzgen_jobs() -> list[BoltzGenJob]:
    """
    Generate BoltzGen jobs with varied binder length ranges.
    Shorter peptides (20-35 aa) resemble semaglutide-class binders.
    Medium (35-60 aa) allows small structured domains.
    Longer (60-90 aa) allows miniproteins.
    """
    jobs: list[BoltzGenJob] = []

    length_configs = [
        ("short_20_30",  20, 30,  "Short peptide range (GLP-1 analog length)"),
        ("short_25_35",  25, 35,  "Short peptide, slightly extended"),
        ("medium_30_45", 30, 45,  "Medium — helix-loop-helix range"),
        ("medium_40_55", 40, 55,  "Medium — small structured domain"),
        ("medium_45_60", 45, 60,  "Medium-long — compact fold"),
        ("long_55_75",   55, 75,  "Long — miniprotein range"),
        ("long_65_85",   65, 85,  "Long miniprotein"),
        ("long_75_90",   75, 90,  "Full miniprotein"),
        # Repeat short range multiple times for sampling diversity
        ("short_20_30_v2", 20, 30, "Short peptide (second sample)"),
        ("short_25_35_v2", 25, 35, "Short peptide extended (second sample)"),
    ]

    for name, lo, hi, desc in length_configs:
        jobs.append(BoltzGenJob(
            name=name,
            binder_length_min=lo,
            binder_length_max=hi,
            description=desc,
            tags=["length_sweep"],
        ))

    return jobs
