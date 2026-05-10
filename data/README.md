# Structure Files

The PDB/CIF files are public data from the RCSB Protein Data Bank and are not stored in this repository.

Download them with:

```bash
mkdir -p Hackathon_Files/Protein_Binder_Hackathon_Files
mkdir -p Hackathon_Files/Small_Molecule_Hackathon_Files

# GLP-1R + Semaglutide complex (Track A target)
curl -o Hackathon_Files/Protein_Binder_Hackathon_Files/7ki0.cif \
  "https://files.rcsb.org/download/7KI0.cif"

curl -o Hackathon_Files/Protein_Binder_Hackathon_Files/7ki0.pdb \
  "https://files.rcsb.org/download/7KI0.pdb"

# GLP-1R + Small molecule (Track B target)
curl -o Hackathon_Files/Small_Molecule_Hackathon_Files/6xox.cif \
  "https://files.rcsb.org/download/6XOX.cif"
```

Or download manually from:
- **7KI0** — https://www.rcsb.org/structure/7KI0 (GLP-1R bound to semaglutide analog + Gs protein)
- **6XOX** — https://www.rcsb.org/structure/6XOX (GLP-1R bound to small molecule agonist LY3502970)

## Chain Key

| Structure | Chain | Contents |
|-----------|-------|----------|
| 7KI0 | F | GLP-1R (entity 6, label_asym_id F) — **use this for BoltzGen** |
| 7KI0 | R | GLP-1R (in pre-extracted `7KI0_GLP1R_chainR.pdb`) — **use this for BindCraft** |

> **Note:** BoltzGen requires the full RCSB mmCIF (`7ki0.cif`) because it validates
> `_entity_poly_seq.entity_id` which is absent from per-chain exports. Pass
> `target_chain_id="F"` when submitting BoltzGen jobs.
