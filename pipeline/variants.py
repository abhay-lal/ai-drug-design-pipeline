"""
Generate peptide variants from the semaglutide baseline sequence.

Strategy:
  1. Single-point mutations  — one AA substitution at each position
  2. C-terminal extensions   — add 1–3 residues to improve stability
  3. N-terminal truncations  — trim leading residues that may be dispensable
  4. Known GLP-1R binding motifs — literature-guided hotspot mutations
  5. Conservative charged substitutions — improve solubility / pI
"""

import itertools
import random
from dataclasses import dataclass, field
from typing import Generator

# ── Constants ────────────────────────────────────────────────────────────────

AMINO_ACIDS = list("ACDEFGHIKLMNPQRSTVWY")

# Semaglutide analog (starting peptide from 7KI0 chain B)
BASELINE = "HAEGTFTSDVSSYLEGQAAKEFIAWLVRGRG"

# Positions (0-based) that are known to be critical for GLP-1R binding
# Refs: He et al. 2023, Zhao et al. 2021
CRITICAL_POSITIONS = {0, 1, 3, 5, 7, 12}  # H, A, G, T, S, Y are key contact residues

# Conservative substitution groups (biochemically similar)
CONSERVATIVE_SUBS: dict[str, list[str]] = {
    "A": ["G", "S", "V"],
    "R": ["K", "H"],
    "N": ["Q", "S"],
    "D": ["E", "N"],
    "C": ["S", "T"],
    "Q": ["N", "E"],
    "E": ["D", "Q"],
    "G": ["A", "S"],
    "H": ["R", "K"],
    "I": ["L", "V", "M"],
    "L": ["I", "V", "M"],
    "K": ["R", "Q"],
    "M": ["L", "I"],
    "F": ["Y", "W", "L"],
    "P": ["A", "G"],
    "S": ["T", "A"],
    "T": ["S", "V"],
    "W": ["F", "Y"],
    "Y": ["F", "H", "W"],
    "V": ["I", "L", "A"],
}

# C-terminal extensions that may improve helicity / receptor contact
CTERM_EXTENSIONS = [
    "K",       # single charged cap
    "RG",      # flexible
    "AAA",     # helix-extending alanines
    "KKKK",    # poly-K solubility tag
    "EIAALEK", # helical extension motif
]

# Known GLP-1R-active modifications (literature-inspired)
HOTSPOT_VARIANTS = [
    # Description, mutation dict  {position: new_aa}
    ("A2G",     {1: "G"}),
    ("T7A",     {6: "A"}),
    ("S8A",     {7: "A"}),
    ("Y13L",    {12: "L"}),
    ("A18C",    {17: "C"}),   # potential disulfide
    ("K21R",    {20: "R"}),
    ("E22D",    {21: "D"}),
    ("I23L",    {22: "L"}),
    ("A24G",    {23: "G"}),
    ("W25F",    {24: "F"}),   # conservative aromatic swap
    ("L26I",    {25: "I"}),
    ("V28I",    {27: "I"}),
    ("R29K",    {28: "K"}),
    ("Aib8",    {7: "A"}),    # alpha-methyl alanine analog
]


# ── Dataclass ────────────────────────────────────────────────────────────────

@dataclass
class PeptideVariant:
    sequence: str
    name: str
    strategy: str
    description: str
    tags: list[str] = field(default_factory=list)

    @property
    def length(self) -> int:
        return len(self.sequence)

    def __str__(self) -> str:
        return f"{self.name}: {self.sequence} ({self.length} aa) [{self.strategy}]"


# ── Generators ───────────────────────────────────────────────────────────────

def single_point_mutations(
    baseline: str = BASELINE,
    skip_critical: bool = False,
    conservative_only: bool = False,
) -> Generator[PeptideVariant, None, None]:
    """Yield all single-point mutants (non-critical positions by default)."""
    for pos, original_aa in enumerate(baseline):
        if skip_critical and pos in CRITICAL_POSITIONS:
            continue
        candidates = CONSERVATIVE_SUBS.get(original_aa, []) if conservative_only else AMINO_ACIDS
        for new_aa in candidates:
            if new_aa == original_aa:
                continue
            mutant = list(baseline)
            mutant[pos] = new_aa
            seq = "".join(mutant)
            yield PeptideVariant(
                sequence=seq,
                name=f"{original_aa}{pos+1}{new_aa}",
                strategy="single_point",
                description=f"Position {pos+1}: {original_aa}→{new_aa}",
                tags=["conservative" if conservative_only else "all_sub"],
            )


def cterm_extensions(
    baseline: str = BASELINE,
    extensions: list[str] | None = None,
) -> Generator[PeptideVariant, None, None]:
    """Yield C-terminal extension variants."""
    for ext in (extensions or CTERM_EXTENSIONS):
        seq = baseline + ext
        yield PeptideVariant(
            sequence=seq,
            name=f"CTerm+{ext}",
            strategy="cterm_extension",
            description=f"C-terminal extension: +{ext} ({len(ext)} aa)",
            tags=["extension"],
        )


def nterm_truncations(
    baseline: str = BASELINE,
    min_length: int = 20,
) -> Generator[PeptideVariant, None, None]:
    """Yield N-terminal truncation variants."""
    for trim in range(1, len(baseline) - min_length + 1):
        seq = baseline[trim:]
        yield PeptideVariant(
            sequence=seq,
            name=f"NTrim{trim}",
            strategy="nterm_truncation",
            description=f"Remove first {trim} residue(s)",
            tags=["truncation"],
        )


def hotspot_variants(
    baseline: str = BASELINE,
) -> Generator[PeptideVariant, None, None]:
    """Yield literature-guided hotspot single or multi-point mutants."""
    for name, mutations in HOTSPOT_VARIANTS:
        mutant = list(baseline)
        for pos, aa in mutations.items():
            if pos < len(mutant):
                mutant[pos] = aa
        yield PeptideVariant(
            sequence="".join(mutant),
            name=name,
            strategy="hotspot",
            description=f"Literature-guided: {name}",
            tags=["literature", "hotspot"],
        )


def double_hotspot_combos(
    baseline: str = BASELINE,
    max_combos: int = 30,
) -> Generator[PeptideVariant, None, None]:
    """Yield random pairwise combinations of hotspot mutations."""
    pairs = list(itertools.combinations(HOTSPOT_VARIANTS, 2))
    random.shuffle(pairs)
    for (n1, m1), (n2, m2) in pairs[:max_combos]:
        # Skip if same position
        if set(m1.keys()) & set(m2.keys()):
            continue
        mutant = list(baseline)
        combined: dict[int, str] = {**m1, **m2}
        for pos, aa in combined.items():
            if pos < len(mutant):
                mutant[pos] = aa
        yield PeptideVariant(
            sequence="".join(mutant),
            name=f"{n1}+{n2}",
            strategy="double_hotspot",
            description=f"Combined mutations: {n1} + {n2}",
            tags=["literature", "combo"],
        )


def generate_all_variants(
    baseline: str = BASELINE,
    conservative_only: bool = True,
    include_extensions: bool = True,
    include_truncations: bool = True,
    include_hotspots: bool = True,
    include_combos: bool = True,
    max_combos: int = 20,
) -> list[PeptideVariant]:
    """
    Collect all variants into a deduplicated list.
    Always includes the baseline as the first entry.
    """
    variants: list[PeptideVariant] = [
        PeptideVariant(
            sequence=baseline,
            name="baseline",
            strategy="baseline",
            description="Semaglutide analog (7KI0 chain B)",
            tags=["baseline"],
        )
    ]

    seen_seqs: set[str] = {baseline}

    def _add(gen):
        for v in gen:
            if v.sequence not in seen_seqs:
                seen_seqs.add(v.sequence)
                variants.append(v)

    _add(single_point_mutations(baseline, conservative_only=conservative_only))
    if include_extensions:
        _add(cterm_extensions(baseline))
    if include_truncations:
        _add(nterm_truncations(baseline))
    if include_hotspots:
        _add(hotspot_variants(baseline))
    if include_combos:
        _add(double_hotspot_combos(baseline, max_combos=max_combos))

    return variants
