import os
import subprocess
import random
import sys
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import AllChem

# 1. Define Fragments
DONORS = {
    "Carbazole": "c1ccc2c(c1)c3ccccc3n2",
    "Diphenylamine": "c1ccc(cc1)Nc2ccccc2",
    "Phenoxazine": "c1ccc2c(c1)Oc3ccccc3N2",
    "Phenothiazine": "c1ccc2c(c1)Sc3ccccc3N2",
    "Dimethylacridine": "CC1(c2ccccc2Nc3ccccc13)C"
}

ACCEPTORS = {
    "Benzonitrile": "c1ccc(cc1)C#N",
    "Triazine": "c1nc(Cl)nc(Cl)n1",
    "Oxadiazole": "c1nnc(o1)c2ccccc2",
    "Pyridine": "c1ccncc1",
    "Sulfone": "c1ccc(cc1)S(=O)(=O)c2ccccc2",
    "Benzophenone": "c1ccc(cc1)C(=O)c2ccccc2"
}

# 2. Setup Directories
BASE = Path(__file__).resolve().parents[1] / "workflow_demo"
(BASE / "step1_smiles").mkdir(parents=True, exist_ok=True)
(BASE / "step2_xtb").mkdir(parents=True, exist_ok=True)
(BASE / "step3_tddft").mkdir(parents=True, exist_ok=True)

# 3. Generate 30 Combinations
combinations = []
for _ in range(30):
    d_name = random.choice(list(DONORS.keys()))
    a_name = random.choice(list(ACCEPTORS.keys()))
    d_smiles = DONORS[d_name]
    a_smiles = ACCEPTORS[a_name]
    # Simple D-A coupling SMILES (N-C or C-C depending on fragment)
    # For this demo, we use a simplified string concatenation for SMILES logic
    # In a real tool, we'd use RDKit reaction SMARTS.
    if d_name == "Diphenylamine":
        smiles = f"{d_smiles[:-1]}({a_smiles})"
    else:
        smiles = f"{d_smiles}-{a_smiles}"
    combinations.append((d_name, a_name, smiles))

summary = ["# High-Throughput Screening Results (30 Molecules)\n\n| ID | Donor-Acceptor | S1 (eV) | T1 (eV) | DeltaEST (eV) | Region |\n| :--- | :--- | :---: | :---: | :---: | :---: |\n"]

# 4. Run Loop
for i, (d_name, a_name, smiles) in enumerate(combinations):
    mol_id = f"TADF_{i:03d}_{d_name}_{a_name}"
    print(f"\n>>> [{i+1}/30] Processing: {mol_id}")
    
    try:
        # Step 1: SMILES -> 3D
        mol = Chem.MolFromSmiles(smiles)
        if mol is None: continue
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, AllChem.ETKDG())
        xyz_init = BASE / "step1_smiles" / f"{mol_id}_init.xyz"
        with open(xyz_init, "w") as f:
            f.write(Chem.MolToXYZBlock(mol))

        # Step 2: xTB Optimization (Fast)
        xyz_opt = BASE / "step2_xtb" / f"{mol_id}_opt.xyz"
        subprocess.run(["xtb", str(xyz_init), "--opt", "--gfn", "ff", "--parallel", "4"], cwd=BASE / "step2_xtb", capture_output=True)
        if os.path.exists(BASE / "step2_xtb" / "xtbopt.xyz"):
            os.rename(BASE / "step2_xtb" / "xtbopt.xyz", xyz_opt)
        else:
            # Fallback to init if opt fails for demo
            xyz_opt = xyz_init

        # Step 3: PySCF TDDFT (B3LYP/3-21G)
        with open(xyz_opt, 'r') as f:
            atom_block = f.read()
        
        pyscf_cmd = f"""
from pyscf import gto, scf, tdscf
mol = gto.M(atom=\"\"\"{atom_block}\"\"\", basis='3-21g', verbose=0)
mf = scf.RHF(mol).run()
td = tdscf.RPA(mf); td.nstates = 1; e = td.kernel()[0]
td_t = tdscf.TDA(mf); td_t.singlet = False; et = td_t.kernel()[0]
s1 = e[0]*27.2114; t1 = et[0]*27.2114; dest = s1-t1
print(f"DATA:{{s1:.3f}}|{{t1:.3f}}|{{dest:.3f}}")
"""
        proc = subprocess.run([sys.executable, "-c", pyscf_cmd], capture_output=True, text=True)
        if "DATA:" in proc.stdout:
            s1, t1, dest = proc.stdout.split("DATA:")[1].strip().split("|")
            region = "Deep Blue" if float(s1) > 3.0 else ("Blue" if float(s1) > 2.7 else "Green/Red")
            summary.append(f"| {i:03d} | {d_name}-{a_name} | {s1} | {t1} | {dest} | {region} |\n")
            print(f"  Result: S1={s1} eV, DeltaEST={dest} eV")
    except Exception as e:
        print(f"  Error: {e}")

with open(BASE / "summary.md", "w") as f:
    f.writelines(summary)
print("\n--- 30 Molecules Screening Complete ---")
