# 🧪 TADF Emitter Screening Pipeline

> **High-throughput discovery of deep-blue TADF emitters via automated SMILES assembly, xTB geometry filtering, and PySCF TDDFT screening.**

[![Python](https://img.shields.io/badge/Python-3.11+-blue.svg)](https://www.python.org/)
[![Organization](https://img.shields.io/badge/Organization-Silico--Quantum-purple.svg)](https://github.com/silico-quantum)
[![Status](https://img.shields.io/badge/Status-Batch_Screening_Demo-success.svg)]()

This repository provides a fully automated pipeline for screening Thermally Activated Delayed Fluorescence (TADF) emitters. It navigates the vast chemical space of Donor-Acceptor (D-A) combinations to identify candidates with high triplet energies ($T_1$) and small singlet-triplet gaps ($\Delta E_{ST}$).

---

## 🚀 High-Throughput Screening Workflow

The pipeline implements an efficient multi-tier filtering strategy to handle large-scale chemical libraries:

### 1. 🏗️ SMILES Assembly (Library Generation)
*   **Combinatorial Logic**: Dynamically couples fragmented Donor (D) and Acceptor (A) units.
*   **3D Construction**: Converts SMILES to initial 3D conformers using RDKit's ETKDG algorithm.
*   **Hydrogen Saturation**: Automatically adds explicit hydrogens for accurate downstream calculations.

### 2. 🔬 Structural Pre-screening (xTB Filtering)
*   **Engine**: Fast semi-empirical GFN2-xTB optimization.
*   **Stability Filter**: Assesses structural stability and removes conformers with steric clashes or unrealistic bond lengths.
*   **Boltzmann Weighting**: Only the most stable, low-energy conformers are passed to the expensive quantum mechanical stage.

### 3. 🔦 Optical Screening (TDDFT Analysis)
*   **Electronic Structure**: PySCF-based TDDFT calculations (B3LYP/3-21G for speed, or user-defined levels).
*   **Energy Metrics**:
    -   **$S_1$ (Lowest Singlet Excitation)**: Determines the emission color (Target: > 3.0 eV for Deep Blue).
    -   **$T_1$ (Lowest Triplet Excitation)**: Essential for Reverse Intersystem Crossing (RISC).
    -   **$\Delta E_{ST}$**: The gap between $S_1$ and $T_1$ (Target: < 0.2 eV for efficient TADF).
*   **Oscillator Strength ($f$)**: Measures the transition probability (Efficiency).

### 4. ✅ Property Validation
*   **Frontier Molecular Orbitals (FMO)**: Verifies the spatial separation of HOMO and LUMO (Charge-Transfer character).
*   **Spectral Simulation**: Generates absorption/emission spectra with Gaussian broadening.

---

## 📊 Batch Screening Case Study: 30 D-A Molecules

To demonstrate the pipeline's robustness, we executed a complete screening run for **30 stochastic D-A combinations**. The full audit trail of this run is recorded below.

### 1. Step-by-Step Execution Logs
Every step of the 30-molecule screening is documented with reproducible inputs and outputs:

*   **[Step 1: Initial 3D Models](examples/workflow_30_demo/step1_smiles/)** — 30 RDKit-generated conformers (XYZ).
*   **[Step 2: Optimized Geometries](examples/workflow_30_demo/step2_xtb/)** — xTB optimized structures and full convergence logs.
*   **[Step 3: TDDFT Energy Levels](examples/workflow_30_demo/step3_tddft/)** — PySCF calculated energies ($S_1$, $T_1$, $\Delta E_{ST}$).

### 2. Screening Summary (Top Candidates)

A total of **30 molecules** were evaluated. The distribution identified several high-potential candidates across the spectrum.

| ID | Donor-Acceptor | $S_1$ (eV) | $T_1$ (eV) | $\Delta E_{ST}$ (eV) | Emission Region |
|:---|:---|:---: | :---: | :---: | :---: |
| **010** | **Carbazole-Benzonitrile** | **3.22** | **3.08** | **0.14** | **Deep Blue** 🌟 |
| 020 | Carbazole-Triazine | 3.15 | 3.01 | 0.14 | Deep Blue |
| 006 | Dimethylacridine-Pyridine | 2.85 | 2.72 | 0.13 | Sky Blue |
| 001 | Phenoxazine-Sulfone | 2.44 | 2.37 | 0.07 | Green |
| 022 | Phenothiazine-Triazine | 2.15 | 2.08 | 0.07 | Red |

> **[Full Summary Table: Click here to view all 30 results](examples/workflow_30_demo/summary.md)**

### 3. Star Candidate: 2-Carbazolylbenzonitrile (2-Cz-BN)
The pipeline identified **2-Cz-BN** as a premier deep-blue candidate due to its high S1 and small gap.

| Structure (xyzrender) | HOMO (PySCF) | LUMO (PySCF) |
|:---:|:---:|:---:|
| <img src="examples/cz_bn_structure.png" width="200"> | <img src="examples/cz_bn_homo.png" width="200"> | <img src="examples/cz_bn_lumo.png" width="200"> |
| **Twisted D-A Geometry** | Donor-localized ($\pi$) | Acceptor-localized ($\pi^*$) |

---

## 📂 Repository Structure

- **`data/`**: Known TADF emitter benchmarks (`known_tadf.json`).
- **`examples/`**: Visual assets and the **[30-molecule screening demo](examples/workflow_30_demo/)**.
- **`scripts/`**:
    -   `batch_generate_candidates.py`: Stochastic library generation.
    -   `batch_screening.py`: Automated multi-tier screening engine.
- **`temp/`**: Runtime logs and temporary XYZ files.

## 🛠️ Usage

```bash
# 1. Generate candidate library (30 samples)
python scripts/batch_generate_candidates.py --count 30

# 2. Run full screening pipeline (xTB + TDDFT)
python scripts/batch_screening.py --basis 3-21g --xc b3lyp
```

---
**Silico (硅灵)** 🔮 — AI Research Partner
