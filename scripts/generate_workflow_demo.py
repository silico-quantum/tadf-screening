import os
import json
import subprocess
from rdkit import Chem
from rdkit.Chem import AllChem

# Define Fragments
DONORS = {
    "Carbazole": "c1ccc2c(c1)c3ccccc3n2",
    "Diphenylamine": "c1ccc(cc1)Nc2ccccc2",
    "Phenoxazine": "c1ccc2c(c1)Oc3ccccc3N2"
}

ACCEPTORS = {
    "Benzonitrile": "c1ccc(cc1)C#N",
    "Triazine": "c1nc(Cl)nc(Cl)n1",
    "Sulfone": "c1ccc(cc1)S(=O)(=O)c2ccccc2"
}

# Combinations to Demo
DEMO_LIST = [
    ("Carbazole", "Benzonitrile", "c1ccc2c(c1)c3ccccc3n2-c4ccc(cc4)C#N"),
    ("Phenoxazine", "Triazine", "c1ccc2c(c1)Oc3ccccc3N2-c4nc(Cl)nc(Cl)n4"),
    ("Diphenylamine", "Sulfone", "c1ccc(cc1)N(c2ccccc2)c3ccc(cc3)S(=O)(=O)c4ccccc4")
]

# Setup Directories
BASE_DIR = "workflow_demo"
os.makedirs(f"{BASE_DIR}/step1_smiles", exist_ok=True)
os.makedirs(f"{BASE_DIR}/step2_xtb", exist_ok=True)
os.makedirs(f"{BASE_DIR}/step3_tddft", exist_ok=True)

summary_lines = ["# Workflow Demo Results\n\n| Molecule | S1 (eV) | T1 (eV) | DeltaEST (eV) | f | Status |\n| :--- | :---: | :---: | :---: | :---: | :---: |\n"]

for d_name, a_name, smiles in DEMO_LIST:
    mol_id = f"{d_name}_{a_name}"
    print(f"\n>>> Processing: {mol_id}")
    
    # --- STEP 1: RDKit 3D Construction ---
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, AllChem.ETKDG())
    initial_xyz = f"{BASE_DIR}/step1_smiles/{mol_id}_initial.xyz"
    with open(initial_xyz, "w") as f:
        f.write(Chem.MolToXYZBlock(mol))
    print(f"  [Step 1] Initial XYZ saved.")

    # --- STEP 2: xTB Optimization ---
    opt_xyz = f"{BASE_DIR}/step2_xtb/{mol_id}_opt.xyz"
    xtb_log = f"{BASE_DIR}/step2_xtb/{mol_id}_xtb.log"
    print(f"  [Step 2] Optimizing with xTB...")
    
    # Run xTB and capture output
    with open(xtb_log, "w") as log_f:
        subprocess.run(["xtb", f"../../{initial_xyz}", "--opt"], 
                       cwd=f"{BASE_DIR}/step2_xtb", 
                       stdout=log_f, stderr=subprocess.STDOUT)
    
    # xTB creates 'xtbopt.xyz' in its cwd
    xtb_output_file = f"{BASE_DIR}/step2_xtb/xtbopt.xyz"
    if os.path.exists(xtb_output_file):
        os.rename(xtb_output_file, opt_xyz)
        print(f"  [Step 2] Optimization complete.")
    else:
        print(f"  [Step 2] WARNING: Optimization might have failed (no xtbopt.xyz).")
        continue

    # --- STEP 3: PySCF TDDFT ---
    print(f"  [Step 3] Running PySCF TDDFT...")
    
    # We pass the content of the XYZ file to gto.M
    with open(opt_xyz, 'r') as f:
        xyz_content = f.read()

    # Use a subprocess to avoid issues with Python state if many runs
    pyscf_script = f"""
from pyscf import gto, scf, tdscf
mol = gto.M(atom=\"\"\"{xyz_content}\"\"\", basis='3-21g', verbose=0)
mf = scf.RHF(mol).run()
# Singlets
td = tdscf.RPA(mf)
td.nstates = 3
e, xy = td.kernel()
osc = td.oscillator_strength()
# Triplets
td_t = tdscf.TDA(mf)
td_t.singlet = False
et, xyt = td_t.kernel()

s1 = e[0] * 27.2114
t1 = et[0] * 27.2114
dest = s1 - t1
f1 = osc[0]
print(f"RESULT|{{s1:.3f}}|{{t1:.3f}}|{{dest:.3f}}|{{f1:.4f}}")
"""
    try:
        proc = subprocess.run(["python3", "-c", pyscf_script], capture_output=True, text=True)
        output = proc.stdout
        if "RESULT|" in output:
            result_str = output.split("RESULT|")[1].strip()
            s1, t1, dest, f1 = result_str.split("|")
            
            region = "Deep Blue" if float(s1) > 3.0 else "Blue/Green"
            summary_lines.append(f"| {mol_id} | {s1} | {t1} | {dest} | {f1} | {region} |\n")
            
            with open(f"{BASE_DIR}/step3_tddft/{mol_id}_pyscf.log", "w") as log_f:
                log_f.write(output)
            print(f"  [Step 3] TDDFT Results: S1={s1}, DeltaEST={dest}")
        else:
            print(f"  [Step 3] FAILED: {proc.stderr}")
    except Exception as e:
        print(f"  [Step 3] EXCEPTION: {e}")

# Save Summary
with open(f"{BASE_DIR}/summary.md", "w") as f:
    f.writelines(summary_lines)

print("\nAll demo steps completed.")
