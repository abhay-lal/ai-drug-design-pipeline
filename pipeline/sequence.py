"""
Extract and manage GLP-1R and peptide sequences.
"""

from pathlib import Path

THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}

# GLP-1R sequence extracted from 7KI0_GLP1R_chainR.pdb (chain R, res 29-423)
GLP1R_SEQUENCE = (
    "TVSLWETVQKWREYRRQCQRSLTEDPPPATDLFCNRTFDEYACWPDGEPGSFVNVSCPWYLPWASSVPQ"
    "GHVYRFCTAEGLWLQKDNSSLPWRDLSECEESPEEQLLFLYIIYTVGYALSFSALVIASAILLGFRHLHC"
    "TRNYIHLNLFASFILRALSVFIKDAALKWMYSTAAQQHQWDGLLSYQDSLSCRLVFLLMQYCVAANYYWLL"
    "VEGVYLYTLLAFSVFSEQWIFRLYVSIGWGVPLLFVVPWGIVKYLYEDEGCWTRNSNMNYWLIIRLPILFA"
    "IGVNFLIFVRVICIVVSKLKANLDIKCRLAKSTLTLIPLLGTHEVIFAFVMDEHARGTLRFIKLFTELSFT"
    "SFQGLMVAILYCFVNNEVQLEFRKSWERWRLE"
)

# Semaglutide analog from 7KI0 (chain B peptide)
SEMAGLUTIDE_SEQUENCE = "HAEGTFTSDVSSYLEGQAAKEFIAWLVRGRG"

# LY3502970 SMILES (from 6XOX, Track B reference)
LY3502970_SMILES = "Cc1cc(NC(=O)c2cnc(N3CCC[C@@H]3C(N)=O)nc2)ccc1F"

# Orthosteric binding pocket residues (1-indexed in chain A when submitted)
POCKET_CONTACTS = [
    ("A", 152),  # Y152
    ("A", 306),  # W306
    ("A", 364),  # E364
    ("A", 380),  # R380
]


def extract_sequence_from_pdb(pdb_path: Path, chain: str = "R") -> str:
    """Extract amino acid sequence from CA atoms of a given chain."""
    residues: dict[int, str] = {}
    with open(pdb_path) as f:
        for line in f:
            if not line.startswith("ATOM"):
                continue
            if line[13:15].strip() != "CA":
                continue
            if line[21].strip() != chain:
                continue
            resname = line[17:20].strip()
            resnum = int(line[22:26].strip())
            if resnum not in residues:
                residues[resnum] = THREE_TO_ONE.get(resname, "X")
    return "".join(residues[k] for k in sorted(residues))


def renumber_pocket_contacts(pdb_path: Path, chain: str = "R") -> list[tuple[str, int]]:
    """
    Map pocket residue numbers from the PDB (which starts at 29) to
    1-based indices as Boltz-2 expects when the chain is provided as a sequence.
    """
    residues: list[int] = []
    seen: set[int] = set()
    with open(pdb_path) as f:
        for line in f:
            if not line.startswith("ATOM"):
                continue
            if line[13:15].strip() != "CA":
                continue
            if line[21].strip() != chain:
                continue
            resnum = int(line[22:26].strip())
            if resnum not in seen:
                seen.add(resnum)
                residues.append(resnum)

    # Build mapping: original PDB residue number → 1-based index
    pdb_to_idx = {pdb_num: idx + 1 for idx, pdb_num in enumerate(sorted(residues))}

    # Pocket residue numbers as they appear in the PDB
    pocket_pdb_nums = [152, 306, 364, 380]
    result = []
    for pdb_num in pocket_pdb_nums:
        if pdb_num in pdb_to_idx:
            result.append(("A", pdb_to_idx[pdb_num]))
    return result
