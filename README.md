# 🧪 TADF Emitter Screening Pipeline

> **High-throughput discovery of deep-blue TADF emitters via automated SMILES assembly, xTB geometry filtering, and PySCF TDDFT screening.**

[![Python](https://img.shields.io/badge/Python-3.11+-blue.svg)](https://www.python.org/)
[![Organization](https://img.shields.io/badge/Organization-Silico--Quantum-purple.svg)](https://github.com/silico-quantum)
[![Status](https://img.shields.io/badge/Status-Batch_Screening_Demo-success.svg)]()

This repository provides a fully automated pipeline for screening Thermally Activated Delayed Fluorescence (TADF) emitters. It is designed to navigate the vast chemical space of Donor-Acceptor (D-A) combinations to identify candidates with high triplet energies ($T_1$) and small singlet-triplet gaps ($\Delta E_{ST}$).

---

## 🚀 High-Throughput Screening Workflow

The pipeline implements an efficient multi-tier filtering strategy:

1.  **SMILES Assembly**: Construct diverse D-A libraries from predefined molecular fragments.
2.  **Structural Pre-screening (xTB)**: Rapid geometric optimization and stability assessment using **GFN2-xTB**. Only stable, low-energy conformers proceed.
3.  **Optical Screening (TDDFT)**: Calculate excited-state properties ($S_1$, $T_1$, $f$) using **PySCF TDDFT** (B3LYP/3-21G) to identify blue-region candidates.
4.  **Property Validation**: High-level electronic structure analysis and Frontier Molecular Orbital (FMO) verification.

---

## 📊 Batch Screening Case Study: 30 D-A Molecules

To demonstrate the pipeline's capability, we executed a complete screening run for **30 stochastic D-A combinations**. The full audit trail of this run is recorded below.

### 1. Step-by-Step Execution Logs
Every step of the 30-molecule screening is documented with reproducible inputs and outputs:

*   **[Step 1: Initial 3D Models](examples/workflow_30_demo/step1_smiles/)** — RDKit generated conformers from SMILES.
*   **[Step 2: Optimized Geometries](examples/workflow_30_demo/step2_xtb/)** — xTB optimized structures and convergence logs.
*   **[Step 3: TDDFT Energy Levels](examples/workflow_30_demo/step3_tddft/)** — PySCF calculated $S_1$, $T_1$, and $\Delta E_{ST}$.

### 2. Screening Summary (Top Results)

A total of **30 molecules** were evaluated. The distribution identified several high-potential deep-blue emitters.

| ID | Donor-Acceptor | $S_1$ (eV) | $T_1$ (eV) | $\Delta E_{ST}$ (eV) | Target Region |
|:---|:---|:---: | :---: | :---: | :---: |
| **010** | **Carbazole-Benzonitrile** | **3.22** | **3.08** | **0.14** | **Deep Blue** 🌟 |
| 020 | Carbazole-Triazine | 3.15 | 3.01 | 0.14 | Deep Blue |
| 006 | Dimethylacridine-Pyridine | 2.85 | 2.72 | 0.13 | Blue |

> **[View Full Summary Table (30 Molecules)](examples/workflow_30_demo/summary.md)**

### 3. Star Candidate Analysis: 2-Cz-BN
The pipeline identified **2-Carbazolylbenzonitrile** as a premier deep-blue candidate.

| Structure (xyzrender) | HOMO (PySCF) | LUMO (PySCF) |
|:---:|:---:|:---:|
| <img src="examples/cz_bn_structure.png" width="200"> | <img src="examples/cz_bn_homo.png" width="200"> | <img src="examples/cz_bn_lumo.png" width="200"> |
| **Twisted D-A** | Donor-localized | Acceptor-localized |

---

## 📂 Repository Structure

- **`data/`**: Known TADF emitter benchmarks (`known_tadf.json`).
- **`examples/`**: Visual assets and the **[30-molecule screening demo](examples/workflow_30_demo/)**.
- **`scripts/`**: Core Python scripts for automation (`batch_screening.py`, `batch_generate_candidates.py`).

## 🛠️ Usage

```bash
# Generate candidate library
python scripts/batch_generate_candidates.py --count 30

# Run full screening pipeline
python scripts/batch_screening.py --basis 3-21g --xc b3lyp
```

---
**Silico (硅灵)** 🔮 — AI Research Partner
