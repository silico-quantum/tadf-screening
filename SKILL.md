# TADF Screening Skill

High-throughput computational screening of Thermally Activated Delayed Fluorescence (TADF) emitters.

## Pipeline

```
SMILES → RDKit 3D → xTB GFN-FF opt → PySCF TDDFT → Filter → Visualize
```

## Requirements

- Python 3.10+
- PySCF >= 2.5
- RDKit (via conda: `conda install -c conda-forge rdkit`)
- xTB 6.x (`brew install xtb`)
- xyzrender (`pip install xyzrender`)

## Key Lessons Learned

### SMILES Source
- **Always use PubChem CID-verified SMILES** — never hand-write SMILES strings
- PubChem API provides canonical SMILES: `https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{CID}/property/IsomericSMILES/JSON`

### Python Environment
- **RDKit**: Use miniconda python3 (conda install)
- **PySCF**: Use system python3 (pip install)

### xTB GFN-FF Optimization
- Fast (~3-5s/molecule) but can have Fortran runtime errors on some molecules
- **Workaround**: Use `--acc 0.5` flag to reduce numerical instability
- Command: `xtb mol.xyz --gfnff --opt --acc 0.5`

### PySCF TDDFT
- **Method**: TDA (Tamm-Dancoff Approximation) is more stable than full TDDFT
- **Basis**: 3-21G for fast screening; def2-SVP for final candidates
- **Oscillator strength**: `td.oscillator_strength()` is a **METHOD** in PySCF 2.12, not a property
  ```python
  f = td.oscillator_strength()  # Correct
  # f = td.oscillator_strength   # Wrong!
  ```

### T₁ Energy Calculation
- **Use ΔSCF (UKS with spin=2)** for triplet ground state energy
- **Do NOT use TDDFT triplet** — less reliable for screening
- ΔE_ST = E(S₁) - E(T₁)

### Large Molecules
- **Skip PySCF for molecules >35 heavy atoms** — memory/time becomes prohibitive
- Use xTB sTDA or semiempirical methods instead

### Molecular Orbital Visualization with xyzrender
- PySCF cubegen output needs modification for xyzrender compatibility
- **Two changes required**:
  1. Negative atom count (e.g., `-37` instead of `37`)
  2. MO index line after atom coordinates
- Use this transformation script:
  ```python
  # Read cube, modify, write
  with open('homo.cube', 'r') as f:
      lines = f.readlines()
  natoms = int(lines[2].split()[0])
  lines[2] = lines[2].replace(f'{natoms}', f'-{natoms}', 1)
  lines[5] = lines[5].rstrip() + '    1\n'  # MO index
  with open('homo_xyzrender.cube', 'w') as f:
      f.writelines(lines)
  ```
- Render: `xyzrender mol.xyz --mo homo_xyzrender.cube --flat-mo --iso 0.04 --png homo.png`

## Filtering Criteria

For TADF candidates:
- Emission wavelength: 450-550 nm (visible light)
- Singlet-triplet gap: ΔE_ST < 0.3 eV (small gap enables RISC)
- Oscillator strength: f > 0.001 (non-zero radiative rate)

## File Structure

```
tadf-screening/
├── README.md
├── requirements.txt
├── SKILL.md
├── data/
│   └── known_tadf.json      # 12 verified TADF molecules
├── examples/
│   ├── screen_12molecules.py
│   └── figures/
│       ├── screening_results.csv
│       ├── mol_*_struct.png
│       ├── mol_*_homo.png
│       ├── mol_*_lumo.png
│       └── mol_*_spectra.png
```

## References

- Tchapet Njafa et al., arXiv:2511.00922v1 — 747-molecule TADF benchmark
- Adachi et al., Nature 2001, 410, 794 — TADF discovery
- Dias et al., Chem. Rev. 2017, 117, 7019 — Organic TADF design principles
