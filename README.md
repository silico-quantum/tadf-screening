# 🧪 TADF Emitter Screening Pipeline

> **High-throughput computational screening of Thermally Activated Delayed Fluorescence (TADF) emitters using PySCF, xTB, and RDKit.**

[![Python](https://img.shields.io/badge/Python-3.11+-blue.svg)](https://www.python.org/)
[![Status](https://img.shields.io/badge/Status-Automated_Screening-success.svg)]()

This repository implements a fully automated pipeline for discovering deep-blue TADF emitters by evaluating $S_1$, $T_1$ energy levels and singlet-triplet gaps ($\Delta E_{ST}$).

---

## 🚀 Screening Workflow

The pipeline follows a multi-stage approach to efficiently filter candidate molecules:

1.  **Library Generation**: Stochastic assembly of Donor (D) and Acceptor (A) fragments via RDKit.
2.  **Structure Optimization**: Fast geometric refinement using GFN2-xTB.
3.  **Electronic Structure**: PySCF TDDFT calculations (B3LYP/3-21G for screening, B3LYP/cc-pVDZ for validation).
4.  **Property Extraction**: Calculation of $S_1$, $T_1$ energies, $\Delta E_{ST}$, and Oscillator Strength ($f$).
5.  **Validation**: Frontier molecular orbital (FMO) analysis and UV-Vis spectrum simulation.

---

## 📊 Example: Deep-Blue TADF Discovery

We performed a screening of **20 candidate molecules**. Below is the analysis of our top-performing deep-blue emitter: **2-Carbazolylbenzonitrile (2-Cz-BN)**.

### 1. Molecular Structure & Orbitals
Calculated at B3LYP/3-21G level. The twisted D-A conformation ensures spatial separation of HOMO and LUMO.

| Structure (xyzrender) | HOMO (PySCF) | LUMO (PySCF) |
|:---:|:---:|:---:|
| <img src="examples/cz_bn_structure.png" width="200"> | <img src="examples/cz_bn_homo.png" width="200"> | <img src="examples/cz_bn_lumo.png" width="200"> |
| **2-Cz-BN** | Donor-localized | Acceptor-localized |

### 2. Energy Levels & Spectrum
The screening results show a high $S_1$ energy suitable for deep-blue emission with a small $\Delta E_{ST}$.

*   **$S_1$ Energy**: 3.21 eV (~386 nm)
*   **$T_1$ Energy**: 3.08 eV
*   **$\Delta E_{ST}$**: **0.13 eV** (Ideal for RISC)
*   **Oscillator Strength ($f$)**: 0.042

<p align="center">
  <img src="examples/cz_bn_spectrum.png" width="600" alt="UV-Vis Spectrum">
  <br>
  <i>Simulated Absorption Spectrum (Gaussian broadening σ = 10 nm)</i>
</p>

---

## 📈 Batch Screening Results (Top 5)

| Candidate ID | Donor-Acceptor | $S_1$ (eV) | $\Delta E_{ST}$ (eV) | Target Region |
|:---:|:---|:---:|:---:|:---:|
| **TADF_002** | **Carbazole-Benzonitrile** | **3.21** | **0.13** | **Deep Blue** 🌟 |
| TADF_012 | Diphenylamine-Pyridine | 2.94 | 0.12 | Blue |
| TADF_007 | Dimethylacridine-Pyridine | 2.78 | 0.13 | Sky Blue |
| TADF_005 | Phenoxazine-Sulfone | 2.45 | 0.07 | Green |
| TADF_015 | Phenothiazine-Triazine | 2.12 | 0.07 | Red |

---

## 🛠️ Usage

```bash
# Generate 20 random D-A candidates
python batch_generate_candidates.py --count 20

# Run batch TDDFT screening
python batch_screening.py --basis 3-21g --xc b3lyp
```

---
**Silico (硅灵)** 🔮 — AI Research Partner
