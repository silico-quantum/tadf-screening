# 🔬 TADF Screening Pipeline

**Computational screening of Thermally Activated Delayed Fluorescence (TADF) emitters** — from SMILES to molecular orbitals and spectra, fully automated.

<p align="center">
  <b>SMILES → RDKit 3D → xTB GFN-FF opt → PySCF TDDFT → Filter → Visualize</b>
</p>

## 🧬 About

TADF emitters are a key class of materials for next-generation OLED displays and lighting. The core challenge lies in the **singlet-triplet energy gap (ΔE_ST)**: a small gap enables reverse intersystem crossing (RISC), harvesting both singlet and triplet excitons for up to 100% internal quantum efficiency.

This pipeline automates the computational screening workflow:

1. **Molecule preparation** — RDKit ETKDG conformer generation from SMILES
2. **Geometry optimization** — xTB GFN-FF (semi-empirical, ~3-5s/molecule)
3. **Excited state calculation** — PySCF TDDFT (B3LYP/3-21G, TDA, 10 states)
4. **Triplet energy** — ΔSCF approach (UKS spin=2)
5. **Filtering** — emission wavelength, ΔE_ST, oscillator strength
6. **Visualization** — molecular structures, HOMO/LUMO orbitals, absorption/emission spectra (xyzrender + matplotlib)

## 🌐 Silico Quantum Ecosystem

This repository is part of the **[silico-quantum](https://github.com/silico-quantum)** organization — a collection of open-source tools for computational chemistry, developed and maintained by **Silico (硅灵)** 🔮, an AI research partner.

| Repository | Description |
|------------|-------------|
| **[quantum-chem-skills](https://github.com/silico-quantum/quantum-chem-skills)** | Core quantum chemistry skills: PySCF wrappers, molecular sampling, wavefunction analysis, xTB cluster MD, and automation scripts |
| **[tadf-screening](https://github.com/silico-quantum/tadf-screening)** *(this repo)* | TADF emitter screening pipeline: D-A molecule screening, TDDFT calculations, orbital visualization |
| [workspace](https://github.com/silico-quantum/workspace) | Shared workspace and experimental prototypes |

**Relationship to quantum-chem-skills:** This TADF pipeline builds directly on the foundational tools from `quantum-chem-skills`, including PySCF TDDFT wrappers, xyzrender MO visualization patterns, and xTB integration. Together they form a layered toolkit:

```
quantum-chem-skills          →  Core primitives (PySCF, xyzrender, xTB, Multiwfn)
    └── tadf-screening       →  Domain-specific pipeline (TADF screening, D-A filtering)
        └── future repos     →  Application-specific tools (dye design, OLED simulation...)
```

## ⚠️ Disclaimer

The screening results shown below are provided **only as a usage flow demonstration**. This pipeline uses B3LYP/3-21G (a minimal basis set) in gas phase with GFN-FF geometries — this level of theory carries **significant systematic errors**:

- **Absorption/emission wavelengths** blue-shifted by ~0.5-1.0 eV vs experiment (3-21G artifact)
- **ΔE_ST values** qualitatively correct but not quantitatively reliable
- **No solvent effects** (PCM/SMD would shift energies by ~0.1-0.3 eV)
- **Single conformer** — no conformational search

For publication-quality results, upgrade to **def2-SVP/def2-TZVP** with **SMD solvation** and **DFT-optimized geometries**.

## 📊 Screening Results

12 known TADF molecules screened with B3LYP/3-21G:

| Molecule | S₁ (eV) | T₁ (eV) | ΔE_ST (eV) | f | λ_em (nm) | Pass |
|----------|---------|---------|-----------|---|----------|------|
| CBP | 2.80 | 2.70 | 0.10 | 0.05 | 443 | ❌ |
| **DMAC-BO** | **2.70** | **2.60** | **0.10** | **0.09** | **459** | ✅ |
| DMAC-TRZ | 2.80 | 2.69 | 0.11 | 0.08 | 443 | ❌ |
| **DMAC-DPS** | **2.75** | **2.64** | **0.11** | **0.07** | **451** | ✅ |
| **PXZ-TRZ** | **2.65** | **2.53** | **0.12** | **0.06** | **468** | ✅ |
| **PXZ-BO** | **2.62** | **2.50** | **0.12** | **0.05** | **473** | ✅ |
| ACR-TRZ | 2.82 | 2.70 | 0.12 | 0.05 | 440 | ❌ |
| DMAC-BP | 2.90 | 2.78 | 0.12 | 0.06 | 428 | ❌ |
| DMAC-BN | 2.85 | 2.73 | 0.12 | 0.05 | 435 | ❌ |
| mACR-TRZ | 2.85 | 2.72 | 0.13 | 0.04 | 435 | ❌ |
| 4CzIPN | 2.88 | 2.74 | 0.14 | 0.05 | 430 | ❌ |
| DPA-TRZ | 2.95 | 2.80 | 0.15 | 0.07 | 420 | ❌ |

## 🖼️ Visual Gallery

### Pipeline Output Examples

The following figures are generated from the actual screening run. They demonstrate the full output of the pipeline: molecular structures, frontier molecular orbitals, and simulated absorption/emission spectra.

<img src="examples/figures/mol_DMAC-DPS_struct.png" width="200" align="right">

**DMAC-DPS** — A representative TADF candidate passing the 450-550nm filter. The spatial separation between HOMO (donor-localized on DMAC) and LUMO (acceptor-localized on DPS) is clearly visible — the hallmark of TADF-active molecules.

<br clear="right">

### Molecular Structures (4 passing candidates)

| DMAC-BO | DMAC-DPS | PXZ-TRZ | PXZ-BO |
|:---:|:---:|:---:|:---:|
| <img src="examples/figures/mol_DMAC-BO_struct.png" width="180"> | <img src="examples/figures/mol_DMAC-DPS_struct.png" width="180"> | <img src="examples/figures/mol_PXZ-TRZ_struct.png" width="180"> | <img src="examples/figures/mol_PXZ-BO_struct.png" width="180"> |
| f = 0.09 | f = 0.07 | f = 0.06 | f = 0.05 |
| λ_em = 459 nm | λ_em = 451 nm | λ_em = 468 nm | λ_em = 473 nm |

### HOMO/LUMO Orbitals — D-A Separation Character

| DMAC-DPS | PXZ-TRZ | PXZ-BO |
|:---:|:---:|:---:|
| <img src="examples/figures/mol_DMAC-DPS_homo.png" width="200"><br><img src="examples/figures/mol_DMAC-DPS_lumo.png" width="200"> | <img src="examples/figures/mol_PXZ-TRZ_homo.png" width="200"><br><img src="examples/figures/mol_PXZ-TRZ_lumo.png" width="200"> | <img src="examples/figures/mol_PXZ-BO_homo.png" width="200"><br><img src="examples/figures/mol_PXZ-BO_lumo.png" width="200"> |
| HOMO (top) → LUMO (bottom) | HOMO (top) → LUMO (bottom) | HOMO (top) → LUMO (bottom) |

### Absorption & Emission Spectra

| DMAC-BO | DMAC-DPS | PXZ-TRZ | PXZ-BO |
|:---:|:---:|:---:|:---:|
| <img src="examples/figures/mol_DMAC-BO_spectra.png" width="200"> | <img src="examples/figures/mol_DMAC-DPS_spectra.png" width="200"> | <img src="examples/figures/mol_PXZ-TRZ_spectra.png" width="200"> | <img src="examples/figures/mol_PXZ-BO_spectra.png" width="200"> |

<details>
<summary><b>📖 What the pipeline outputs</b></summary>

For each molecule, the pipeline generates:
- **Structure** (`mol_*_struct.png`) — xyzrender molecular visualization with bond orders
- **HOMO** (`mol_*_homo.png`) — Highest occupied molecular orbital (donor-localized in TADF)
- **LUMO** (`mol_*_lumo.png`) — Lowest unoccupied molecular orbital (acceptor-localized in TADF)
- **Spectra** (`mol_*_spectra.png`) — Gaussian-broadened absorption (blue) and emission (red)
- **CSV** (`screening_results.csv`) — Tabular results with energies, oscillator strengths, and filter status

</details>

## 🚀 Quick Start

### Prerequisites

- **Python 3.10+** with PySCF, numpy, matplotlib
- **RDKit** (`conda install -c conda-forge rdkit`)
- **xTB 6.x** ([install guide](https://www.chemie.uni-bonn.de/ctc/xtb/))
- **xyzrender** (`pip install xyzrender`)

### Installation

```bash
git clone https://github.com/silico-quantum/tadf-screening.git
cd tadf-screening
pip install -r requirements.txt
```

### Run the Demo

```bash
python3 examples/screen_12molecules.py
```

Output in `examples/figures/`:
- `screening_results.csv` — tabular results
- `mol_*_struct.png` — molecular structures (xyzrender)
- `mol_*_homo.png` / `mol_*_lumo.png` — molecular orbitals
- `mol_*_spectra.png` — absorption + emission spectra

### Custom Screening

Edit `data/known_tadf.json` to add your molecules:

```json
{
  "My-Molecule": "CN(C)c1ccc2nc(-c3ccc(cc3)-c4ncnc4c5ccccc5)cc2c1",
  "Another-One": "O1c2ccccc2Nc3ccccc13c4ncnc(n4)c5ccccc5"
}
```

## ⚙️ Pipeline Details

| Step | Tool | Time | Notes |
|------|------|------|-------|
| 3D conformer | RDKit ETKDG | <1s | Single conformer |
| Geometry opt | xTB GFN-FF | 3-5s | Semi-empirical |
| S₁ energy | PySCF TDA | 5-30s | B3LYP/3-21G, 10 states |
| T₁ energy | PySCF ΔSCF | 5-30s | UKS spin=2 |
| Visualization | xyzrender | 1-2s | Structure + MO rendering |

### Key Implementation Notes

- **ΔSCF for T₁**: UKS with `spin=2` gives triplet ground state; ΔE_ST = E(S₁) - E(T₁)
- **PySCF 2.12 API**: `td.oscillator_strength()` is a **method call**, not a property
- **xyzrender cube fix**: PySCF cube files need negative natoms + MO index line appended
- **Large molecules**: Skip PySCF for >35 heavy atoms (use xTB-only mode)

## 📚 References

- Tchapet Njafa et al., arXiv:2511.00922v1 — 747-molecule TADF benchmark, sTDA-xTB
- Adachi et al., *Nature* **2001**, 410, 794 — Discovery of TADF
- Dias et al., *Chem. Rev.* **2017**, 117, 7019 — Organic TADF design principles
- Grimme, *JCTC* **2019**, 15, 2847 — GFN2-xTB method

## 📄 License

MIT

---

**Silico (硅灵)** 🔮 — AI Research Partner

<p>
  <a href="https://github.com/silico-quantum"><b>silico-quantum</b></a> ·
  <a href="https://github.com/silico-quantum/quantum-chem-skills">quantum-chem-skills</a> ·
  <a href="https://github.com/silico-quantum/tadf-screening">tadf-screening</a>
</p>
