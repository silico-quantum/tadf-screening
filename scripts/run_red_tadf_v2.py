import os
import subprocess
import sys
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import AllChem

# Target: Red-emitting TADF molecule
# Combination: Phenothiazine-Benzophenone (PTZ-BP) - Known for red-shifted emission
mol_id = "Red_TADF_PTZ_BP"
smiles = "c1ccc(cc1)C(=O)c2ccc(cc2)N3c4ccccc4Sc5ccccc35"

BASE = Path(__file__).resolve().parents[1] / "red_tadf_demo"
BASE.mkdir(parents=True, exist_ok=True)

print(f">>> Running Red-TADF workflow for: {mol_id}")

# 1. SMILES -> 3D (RDKit)
mol = Chem.MolFromSmiles(smiles)
mol = Chem.AddHs(mol)
AllChem.EmbedMolecule(mol, AllChem.ETKDG())
AllChem.MMFFOptimizeMolecule(mol)
xyz_init = BASE / f"{mol_id}_init.xyz"

# Clean XYZ extraction for PySCF compatibility
xyz_block = Chem.MolToXYZBlock(mol)
with open(xyz_init, "w") as f:
    f.write(xyz_block)
print(f"  [1/3] Step 1: Initial XYZ saved.")

# 2. xTB Optimization (Fast)
xyz_opt = BASE / f"{mol_id}_opt.xyz"
print(f"  [2/3] Step 2: Optimizing with xTB...")
subprocess.run(["xtb", xyz_init, "--opt", "--parallel", "4"], cwd=BASE, capture_output=True)
if os.path.exists(f"{BASE}/xtbopt.xyz"):
    os.rename(f"{BASE}/xtbopt.xyz", xyz_opt)
    print(f"  [2/3] Step 2: Optimized XYZ saved.")
else:
    xyz_opt = xyz_init
    print(f"  [2/3] Step 2: Fallback to initial XYZ.")

# 3. PySCF TDDFT
print(f"  [3/3] Step 3: Running PySCF TDDFT (B3LYP/3-21G)...")
# Extract atoms for PySCF
lines = open(xyz_opt).readlines()[2:]
atom_pyscf = "".join(lines)

pyscf_script = f"""
from pyscf import gto, scf, tdscf
mol = gto.M(atom=\"\"\"{atom_pyscf}\"\"\", basis='3-21g', verbose=0)
mf = scf.RHF(mol).run()
td = tdscf.RPA(mf); td.nstates = 5; e = td.kernel()[0]
osc = td.oscillator_strength()
td_t = tdscf.TDA(mf); td_t.singlet = False; et = td_t.kernel()[0]
s1 = e[0]*27.2114; t1 = et[0]*27.2114; dest = s1-t1; f1 = osc[0]
print(f"DATA:{{s1:.3f}}|{{t1:.3f}}|{{dest:.3f}}|{{f1:.4f}}")
"""
proc = subprocess.run([sys.executable, "-c", pyscf_script], capture_output=True, text=True)
if "DATA:" in proc.stdout:
    s1, t1, dest, f1 = proc.stdout.split("DATA:")[1].strip().split("|")
    nm = 1240.0 / float(s1)
    print(f"  [3/3] Result: S1={s1} eV ({nm:.1f} nm), DeltaEST={dest} eV, f={f1}")
    with open(BASE / "results.md", "w") as f:
        f.write(f"# Red TADF Result\n- Molecule: {mol_id}\n- S1: {s1} eV ({nm:.1f} nm)\n- T1: {t1} eV\n- DeltaEST: {dest} eV\n- f: {f1}\n")
else:
    print(f"  [3/3] FAILED: {proc.stderr}")
