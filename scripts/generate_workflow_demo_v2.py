import os
import subprocess
import sys
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import AllChem

# Define Fragments (D & A)
DEMO_LIST = [
    ("Carbazole_Benzonitrile", "c1ccc2c(c1)c3ccccc3n2-c4ccc(cc4)C#N"),
    ("Phenoxazine_Triazine", "c1ccc2c(c1)Oc3ccccc3N2-c4nc(Cl)nc(Cl)n4"),
    ("Diphenylamine_Sulfone", "c1ccc(cc1)N(c2ccccc2)c3ccc(cc3)S(=O)(=O)c4ccccc4")
]

# Workdir setup
BASE = Path(__file__).resolve().parents[1] / "workflow_demo"
(BASE / "step1_smiles").mkdir(parents=True, exist_ok=True)
(BASE / "step2_xtb").mkdir(parents=True, exist_ok=True)
(BASE / "step3_tddft").mkdir(parents=True, exist_ok=True)

summary = ["# Workflow Demo Results\n\n| Molecule | S1 (eV) | T1 (eV) | DeltaEST (eV) | f | Region |\n| :--- | :---: | :---: | :---: | :---: | :---: |\n"]

for mol_id, smiles in DEMO_LIST:
    print(f"\n>>> Running workflow for: {mol_id}")
    
    # 1. SMILES -> 3D (RDKit)
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, AllChem.ETKDG())
    xyz_init = BASE / "step1_smiles" / f"{mol_id}_init.xyz"
    with open(xyz_init, "w") as f:
        f.write(Chem.MolToXYZBlock(mol))
    print(f"  [1/3] Step 1: Initial XYZ saved to {xyz_init}")

    # 2. xTB Optimization
    xyz_opt = BASE / "step2_xtb" / f"{mol_id}_opt.xyz"
    print(f"  [2/3] Step 2: Optimizing with xTB...")
    # Run in step2 dir to keep it clean
    subprocess.run(["xtb", str(xyz_init), "--opt", "--parallel", "4"], cwd=BASE / "step2_xtb", capture_output=True)
    if os.path.exists(BASE / "step2_xtb" / "xtbopt.xyz"):
        os.rename(BASE / "step2_xtb" / "xtbopt.xyz", xyz_opt)
        print(f"  [2/3] Step 2: Optimized XYZ saved to {xyz_opt}")
    else:
        print(f"  [2/3] Step 2: FAILED (xTB optimization did not produce xtbopt.xyz)")
        continue

    # 3. PySCF TDDFT
    print(f"  [3/3] Step 3: Running PySCF TDDFT (B3LYP/3-21G)...")
    with open(xyz_opt, 'r') as f:
        atom_block = f.read()
    
    pyscf_cmd = f"""
from pyscf import gto, scf, tdscf
import numpy as np
mol = gto.M(atom=\"\"\"{atom_block}\"\"\", basis='3-21g', verbose=0)
mf = scf.RHF(mol).run()
td = tdscf.RPA(mf)
td.nstates = 3
e, xy = td.kernel()
osc = td.oscillator_strength()
td_t = tdscf.TDA(mf)
td_t.singlet = False
et, xyt = td_t.kernel()
s1 = e[0]*27.2114; t1 = et[0]*27.2114; dest = s1-t1; f1 = osc[0]
print(f"DATA:{{s1:.3f}}|{{t1:.3f}}|{{dest:.3f}}|{{f1:.4f}}")
"""
    proc = subprocess.run([sys.executable, "-c", pyscf_cmd], capture_output=True, text=True)
    if "DATA:" in proc.stdout:
        res = proc.stdout.split("DATA:")[1].strip().split("|")
        s1, t1, dest, f1 = res
        region = "Deep Blue" if float(s1) > 3.0 else ("Blue" if float(s1) > 2.7 else "Green/Yellow")
        summary.append(f"| {mol_id} | {s1} | {t1} | {dest} | {f1} | {region} |\n")
        print(f"  [3/3] Step 3: SUCCESS. S1={s1} eV, DeltaEST={dest} eV")
    else:
        print(f"  [3/3] Step 3: FAILED. Error: {proc.stderr}")

with open(BASE / "summary.md", "w") as f:
    f.writelines(summary)

print("\n--- Workflow Demo Complete ---")
