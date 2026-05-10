"""
Build Boltz-2 YAML job files for Track A peptide binder submissions.
"""

from pathlib import Path
from typing import Any

import yaml

from sequence import GLP1R_SEQUENCE, POCKET_CONTACTS
from variants import PeptideVariant

# Remapped pocket contacts (1-based indices matching the sequence provided)
# Y152→117, W306→271, E364→325, R380→341  (PDB res 29 = index 1)
REMAPPED_POCKET_CONTACTS: list[tuple[str, int]] = [
    ("A", 117),
    ("A", 271),
    ("A", 325),
    ("A", 341),
]


def build_boltz2_yaml(
    variant: PeptideVariant,
    receptor_sequence: str = GLP1R_SEQUENCE,
    pocket_contacts: list[tuple[str, int]] = REMAPPED_POCKET_CONTACTS,
) -> dict[str, Any]:
    """Return a Boltz-2 YAML dict for a peptide binder job."""
    return {
        "version": 1,
        "sequences": [
            {
                "protein": {
                    "id": "A",
                    "sequence": receptor_sequence,
                }
            },
            {
                "protein": {
                    "id": "B",
                    "sequence": variant.sequence,
                }
            },
        ],
        "constraints": [
            {
                "pocket": {
                    "binder": "B",
                    "contacts": [[chain, resnum] for chain, resnum in pocket_contacts],
                }
            }
        ],
    }


def write_yaml_job(
    variant: PeptideVariant,
    output_dir: Path,
    receptor_sequence: str = GLP1R_SEQUENCE,
    pocket_contacts: list[tuple[str, int]] = REMAPPED_POCKET_CONTACTS,
) -> Path:
    """Write a YAML job file for the given variant. Returns the file path."""
    output_dir.mkdir(parents=True, exist_ok=True)
    job = build_boltz2_yaml(variant, receptor_sequence, pocket_contacts)
    # Safe filename: replace + and spaces
    safe_name = variant.name.replace("+", "_").replace(" ", "_")
    path = output_dir / f"{safe_name}.yaml"
    with open(path, "w") as f:
        yaml.dump(job, f, default_flow_style=False, sort_keys=False, allow_unicode=True)
    return path


def write_all_yaml_jobs(
    variants: list[PeptideVariant],
    output_dir: Path,
    receptor_sequence: str = GLP1R_SEQUENCE,
) -> list[Path]:
    """Write YAML job files for all variants. Returns list of paths."""
    paths = []
    for v in variants:
        p = write_yaml_job(v, output_dir, receptor_sequence)
        paths.append(p)
    return paths
