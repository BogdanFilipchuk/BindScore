<div align="center">
  <img src="Notebooks/logo.jpg" width="220"/>
</div>

# BindScore — Protein-Protein Binding Free Energy Estimator

> A static-structure, physics-motivated scorer that estimates ΔH, ΔS, and ΔG for protein-protein complexes directly from a PDB file — no molecular dynamics engine required.

---

## Authors

| Name | Email | Affiliation |
|------|-------|-------------|
| Bogdan Filipchuk | bogdan.filipchuk@epfl.ch | EPFL |
| Tuna Karasu | tuna.karasu@epfl.ch | EPFL |
| Andrey Babenko | andrey.babenko@epfl.ch | EPFL |
| Maximus van den Bogaard | maximus.vandenbogaard@epfl.ch | EPFL |

---

## Overview

Predicting protein–protein binding affinity from structure is a long-standing problem in computational biology. It has been addressed by a wide spectrum of tools — from MD-based methods (FEP, MM-PBSA) through empirical scorers (FoldX, PRODIGY, Rosetta InterfaceAnalyzer) to modern deep-learning predictors (AlphaFold-Multimer-derived models, DeepDDG). Each family trades accuracy against speed, transparency, and the amount of structural input required.

BindScore is a from-scratch implementation in the spirit of fast empirical scorers, written as a course project to make every scoring term physically explicit rather than to compete with production tools. The baseline design goal was a running time measured in seconds — no MD engine is used.

The tool scores two entities (protein chains, small molecules, or metals) inside a single PDB file and returns a per-component breakdown of ΔH and ΔS, from which an estimate of ΔG can be assembled.

---

## Scientific Background

From the beginning, a clear separation was made between the enthalpic and entropic sides of the project.

### Enthalpy (ΔH)

The enthalpic pipeline detects atom-atom contacts at the interface, classifies each one into a physical interaction type, and scores it with an analytic potential that fits that type:

| Interaction type | Potential |
|-----------------|-----------|
| Hydrogen bond | Gaussian well centered at 3.0 Å, depth 8 kJ/mol |
| Hydrophobic contact | Gaussian well centered at 3.8 Å, depth 0.5 kJ/mol |
| Salt bridge | Screened Coulomb (ε_r = 50), aggregated per residue pair |
| π–π stacking | LJ-like well, depth and eq. distance depend on ring-normal angle |
| Halogen / disulfide / metal coordination | LJ wells with type-specific depths |
| Dipole–dipole | Fraction of screened Coulomb between atomic partial charges |
| Desolvation penalty | One water H-bond enthalpy (4.25 kJ/mol) per buried polar contact |

The enthalpic pipeline is complete and produces a per-interaction breakdown.

### Entropy (ΔS) — preliminary

Binding entropy is broken down into four independently computed components:

- **Trans+Rot:** Analytical estimate from the Sackur–Tetrode and rigid-rotor partition functions (Finkelstein–Janin cage-volume convention). Costs microseconds and reproduces the canonical ~6 kcal/mol rigid-body penalty.
- **Solvent (hydrophobic):** Sun (2022) water-partitioning model, which separates interfacial water molecules perturbed by the protein surface from bulk tetrahedral water. The bound complex is compared against the summed free-chain entropies to extract the net solvent reorganization.
- **Sidechain conformational:** Pickett–Sternberg (1993) empirical scale — a per-residue ΔS value for interface residues that lose SASA on binding.
- **Backbone (NMA):** Tidor–Karplus mode-matching approach using ProDy ANM normal modes, with Hungarian assignment between complex and isolated-chain modes.

> **Note:** The entropy modules are preliminary and under active validation. An earlier RR-FEP approach using OpenMM was abandoned because a single 3-atom molecule required tens of minutes of simulation time.

---

## Benchmarks

Enthalpy scoring was validated against 22 complexes with experimental ΔG values, giving a **68.18 % success rate** within a ±20 kJ/mol threshold.

| Reference | PDB ID | Experimental (kJ/mol) | Calculated (kJ/mol) | Error (kJ/mol) | Status |
|:----------|:------:|:---------------------:|:-------------------:|:--------------:|:------:|
| Xie et al | **1DPU** | -70.32 | -24.42 | +45.90 | FAIL |
| Schmidt et al | **1RST** | -52.58 | -36.27 | +16.31 | Pass |
| Liu and Vogel | **2LQC** | -28.92 | -27.11 | +1.81 | Pass |
| Yu et al | **2MNU** | -19.25 | -40.65 | -21.40 | FAIL |
| Grace et al | **2MWY** | -73.25 | -60.14 | +13.11 | Pass |
| Eulitz et al | **4F14** | -35.41 | -75.01 | -39.60 | FAIL |
| Tallant et al | **4Q6F** | -41.06 | -43.15 | -2.09 | Pass |
| Clark et al | **5E0M** | -34.32 | -36.29 | -1.97 | Pass |
| Ponna et al | **5OVC** | -28.46 | -37.06 | -8.60 | Pass |
| Murthy et al | **6EVO** | -36.42 | -46.81 | -10.39 | Pass |
| Huber et al | **6H8C** | -24.74 | -65.29 | -40.55 | FAIL |
| Frisch, Schreiber, Johnson | **1BRS** | -77.85 | -61.24 | +16.61 | Pass |
| Filippakopoulos et al | **3MXF** | -35.24 | -58.82 | -23.58 | FAIL |
| Filippakopoulos et al | **3U5L** | -25.78 | -44.96 | -19.18 | Pass |
| Lucas et al | **4LZR** | -37.67 | -42.25 | -4.58 | Pass |
| Gacias et al | **4QB3** | -27.71 | -50.82 | -23.11 | FAIL |
| Picaud et al | **4XY9** | -25.49 | -36.19 | -10.70 | Pass |
| Xue et al | **5DOC** | -42.70 | -40.00 | +2.70 | Pass |
| Hügle et al | **5D3S** | -40.90 | -41.72 | -0.82 | Pass |
| Raux et al | **5DW2** | -42.28 | -44.52 | -2.24 | Pass |
| Montenegro et al | **5FBX** | -65.18 | -46.13 | +19.05 | Pass |
| Picaud et al | **5IGK** | -46.42 | -55.27 | -8.85 | Pass |

Failures are predominantly caused by the static-structure assumption: the scorer operates on a single frozen conformation and cannot account for conformational selection or induced fit.

---

## Installation

```bash
git clone https://github.com/BogdanFilipchuk/BindScore.git
cd "BindScore Project"
pip install -e .
```

**Requirements:** Python ≥ 3.10. All dependencies (biopython, numpy, scipy, prody, freesasa, click, streamlit, …) are declared in `pyproject.toml` and installed automatically.

---

## CLI Usage

After installation the `bindscore` command is available in your environment.

### Inspect a PDB structure

```bash
bindscore getchains 1BRS
```
Fetches the structure from RCSB if not local and prints all polypeptide chains and small molecules.

### List atom-pair contacts

```bash
bindscore interactions 1BRS A D
bindscore interactions 1BRS A D --threshold 4.5
```
Classifies and prints every detected contact between two chains.

### Score a complex

```bash
bindscore score 1BRS A D
bindscore score 1BRS A D --breakdown
bindscore score 1BRS A D --breakdown --json
```

| Option | Short | Description |
|--------|-------|-------------|
| `--temperature` | `-T` | Temperature in Kelvin (default 300) |
| `--breakdown` | `-b` | Print per-component ΔS breakdown |
| `--json` | | Output results as JSON |
| `--show_failed` | `-sf` | Show detailed errors from failed submodules |

### Download a PDB file

```bash
bindscore fetch 1BRS --out ./pdbs/
```

### Launch the web interface

```bash
bindscore app
bindscore app --port 8080
```

---

## Python API

```python
from bindscore.pdb_file_treatment.pdb_utils_fetch import fetch_pdb_data
from bindscore.parsing.pdb_utils_protein import Protein_Structure
from bindscore.parsing.pdb_utils_inter import Interaction
from bindscore.scoring.pdb_utils_enthalpy import binding_enthalpy
from bindscore.scoring import compute_total_entropy

pdb_data = fetch_pdb_data("1BRS")
protein  = Protein_Structure(pdb_data)
inter    = Interaction(protein, "A", "D", threshold=5.0)

dH = binding_enthalpy(inter)["TOTAL"]          # kJ/mol
dS = compute_total_entropy("1BRS.pdb", "A", "D", T=298.15, return_breakdown=True)
```

---

## Project Structure

```
src/bindscore/
├── pdb_file_treatment/   # PDB fetch, parsing, radius estimation
├── parsing/              # Protein_Structure, Interaction, enthalpy helpers
├── scoring/              # ΔS submodules (trans_rot, sidechain, backbone, solvent)
│   └── total_entropy.py  # Orchestrates all four ΔS components
└── web/                  # Streamlit app
```

---

## Limitations

- **Static structure:** one conformation is scored; conformational selection and induced fit are not captured.
- **No protonation prediction:** ionization states of His, Asp, Glu are assumed at pH 7.
- **Entropy is preliminary:** all four ΔS modules pass basic sanity checks but have not been validated against a benchmark set.
- **Water-mediated contacts are not scored:** crystallographic bridging waters are ignored.

---

## AI Policy

During the creation of this project, AI tools (Claude, ChatGPT, etc.) were used respecting the course guidelines. The final version functioning version of the package* is conceived by the authors, with comprehensive research and understanding behind key functionalities. AI tools were used for assistance in research and code conceptualisation.
*By 25/05/2026, the Entropy submodules that are tagged MODULE NOT READY are just sketches. They are to be implemented properly in the future with compliance to the AI policy stated before.

---

## References

- Pickett & Sternberg (1993). *J. Mol. Biol.* 231:825–839. Sidechain rotameric entropy scale.
- Finkelstein & Janin (1989). *Protein Eng.* 3(1):1–3. Trans+rot entropy cage convention.
- Sun (2022). Hydrophobic solvation water-partitioning model.
- Tidor & Karplus (1994). *J. Mol. Biol.* 238:405–414. NMA-based backbone entropy.
- Schreiber & Fersht (1993). *Biochemistry* 32(19):5145–5150. Barnase–Barstar benchmark.
- Nidhi Singh & Arieh Warshel (2010). *Proteins* 78(7):1705–1723. doi:10.1002/prot.22689. Entropy decomposition reference.
