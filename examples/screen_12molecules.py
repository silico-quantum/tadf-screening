#!/usr/bin/env python3
"""
TADF Screening Pipeline - Robust Version
Runs complete screening on 12 known TADF molecules
"""

import json
import os
import subprocess
import sys
import numpy as np
from pathlib import Path
import shutil

# Force unbuffered output
sys.stdout.reconfigure(line_buffering=True)
sys.stderr.reconfigure(line_buffering=True)

# Paths
WORK_DIR = Path(__file__).resolve().parents[1]
DEFAULT_PYTHON = os.environ.get("TADF_SCREENING_PYTHON", sys.executable)
RDKIT_PYTHON = os.environ.get("RDKIT_PYTHON", DEFAULT_PYTHON)
PYSCF_PYTHON = os.environ.get("PYSCF_PYTHON", DEFAULT_PYTHON)
XTB_BIN = os.environ.get("XTB_BIN", shutil.which("xtb") or "xtb")
XYZRENDER_BIN = os.environ.get("XYZRENDER_BIN", shutil.which("xyzrender") or "xyzrender")

# Directories
DATA_FILE = WORK_DIR / "data" / "known_tadf.json"
OUTPUT_DIR = WORK_DIR / "examples" / "figures"
TEMP_DIR = WORK_DIR / "temp"

# Create directories
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
TEMP_DIR.mkdir(parents=True, exist_ok=True)

print("=" * 60, flush=True)
print("TADF Screening Pipeline", flush=True)
print("=" * 60, flush=True)

def count_heavy_atoms(smiles):
    """Count non-hydrogen atoms using RDKit"""
    script = f'''
import sys
from rdkit import Chem
mol = Chem.MolFromSmiles("{smiles}")
if mol is None:
    print(-1)
    sys.exit(1)
mol = Chem.AddHs(mol)
heavy_atoms = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() != 1)
print(heavy_atoms)
'''
    result = subprocess.run([RDKIT_PYTHON, "-c", script], capture_output=True, text=True, timeout=30)
    return int(result.stdout.strip())

def generate_3d_conformer(smiles, name):
    """Generate 3D conformer and write XYZ file"""
    xyz_file = TEMP_DIR / f"{name}_initial.xyz"
    
    script = f'''
import sys
from rdkit import Chem
from rdkit.Chem import AllChem

smiles = "{smiles}"
mol = Chem.MolFromSmiles(smiles)
if mol is None:
    print("ERROR: Invalid SMILES", file=sys.stderr)
    sys.exit(1)

mol = Chem.AddHs(mol)
res = AllChem.EmbedMolecule(mol, AllChem.ETKDG())
if res == -1:
    print("ERROR: Embedding failed", file=sys.stderr)
    sys.exit(1)

AllChem.MMFFOptimizeMolecule(mol)

# Write XYZ
conf = mol.GetConformer()
atoms = mol.GetAtoms()
n_atoms = len(atoms)

with open("{xyz_file}", "w") as f:
    f.write(f"{{n_atoms}}\\n")
    f.write(f"Generated from SMILES: {smiles}\\n")
    for i, atom in enumerate(atoms):
        pos = conf.GetAtomPosition(i)
        symbol = atom.GetSymbol()
        f.write(f"{{symbol}} {{pos.x:.6f}} {{pos.y:.6f}} {{pos.z:.6f}}\\n")

print("SUCCESS")
'''
    result = subprocess.run([RDKIT_PYTHON, "-c", script], capture_output=True, text=True, timeout=60)
    if "SUCCESS" not in result.stdout:
        print(f"  [ERROR] 3D conformer generation failed for {name}: {result.stderr}", flush=True)
        return None
    
    return xyz_file

def run_xtb_optimization(xyz_file, name):
    """Run xTB geometry optimization"""
    opt_xyz = TEMP_DIR / f"{name}_xtbopt.xyz"
    log_file = TEMP_DIR / f"{name}_xtb.log"
    
    # Run xTB with simpler output
    cmd = f"cd {TEMP_DIR} && {XTB_BIN} {xyz_file.name} --gfnff --opt > {log_file.name} 2>&1"
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True, timeout=300)
    
    # Parse optimized geometry from xtbopt.xyz
    xtbopt_file = TEMP_DIR / "xtbopt.xyz"
    if xtbopt_file.exists():
        shutil.move(str(xtbopt_file), str(opt_xyz))
        return opt_xyz
    else:
        print(f"  [WARNING] xTB optimization failed for {name}, using initial geometry", flush=True)
        return xyz_file

def run_pyscf_tddft(xyz_file, name):
    """Run PySCF TDDFT calculation"""
    # Parse XYZ to PySCF format
    with open(xyz_file) as f:
        lines = f.readlines()[2:]  # Skip header
    
    atoms = []
    for line in lines:
        parts = line.split()
        if len(parts) >= 4:
            symbol, x, y, z = parts[0], float(parts[1]), float(parts[2]), float(parts[3])
            atoms.append(f"{symbol} {x:.6f} {y:.6f} {z:.6f}")
    
    pyscf_atoms = "; ".join(atoms)
    
    # Create PySCF script
    script = f'''
import sys
import json
from pyscf import gto, dft, tdscf
import numpy as np

try:
    # Build molecule
    mol = gto.M(
        atom="{pyscf_atoms}",
        basis="3-21g",
        charge=0,
        spin=0,
        symmetry=False,
        verbose=0
    )
    
    # Ground state DFT
    mf = dft.RKS(mol, xc="b3lyp")
    mf.kernel()
    e_s0 = mf.e_tot
    
    # TDDFT for excited states (TDA)
    td = tdscf.TDA(mf)
    td.nstates = 10
    td.kernel()
    
    # S1 energy and oscillator strength
    e_s1_hartree = td.e[0]
    e_s1_ev = e_s1_hartree * 27.2114
    
    # Oscillator strength (method call in PySCF 2.12)
    osc_strengths = td.oscillator_strength()
    osc = osc_strengths[0] if len(osc_strengths) > 0 else 0.0
    
    # Triplet T1 via ΔSCF
    mol_t = gto.M(
        atom="{pyscf_atoms}",
        basis="3-21g",
        charge=0,
        spin=2,
        symmetry=False,
        verbose=0
    )
    mf_t = dft.UKS(mol_t, xc="b3lyp")
    mf_t.kernel()
    e_t1_hartree = mf_t.e_tot
    
    e_t1_ev = (e_t1_hartree - e_s0) * 27.2114
    delta_est = e_s1_ev - e_t1_ev
    
    result = {{
        "e_s1_ev": float(e_s1_ev),
        "e_t1_ev": float(e_t1_ev),
        "delta_est_ev": float(delta_est),
        "osc": float(osc),
        "success": True
    }}
    
    print(json.dumps(result))
    
except Exception as e:
    result = {{"error": str(e), "success": False}}
    print(json.dumps(result))
    sys.exit(1)
'''
    
    result_file = TEMP_DIR / f"{name}_pyscf_result.json"
    cmd = f"{PYSCF_PYTHON} -c '{script}'"
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True, timeout=600)
    
    try:
        data = json.loads(result.stdout.strip())
        if data.get("success"):
            return data["e_s1_ev"], data["e_t1_ev"], data["osc"]
        else:
            print(f"  [ERROR] PySCF failed for {name}: {data.get('error', 'Unknown error')}", flush=True)
            return None, None, None
    except Exception as e:
        print(f"  [ERROR] Failed to parse PySCF output for {name}: {e}", flush=True)
        print(f"  Output: {result.stdout[:200]}", flush=True)
        return None, None, None

def generate_structure_visualization(xyz_file, name):
    """Generate molecular structure visualization"""
    output_file = OUTPUT_DIR / f"mol_{name}_struct.png"
    cmd = f"{XYZRENDER_BIN} {xyz_file} -o {output_file} -t --bo -S 600"
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    
    if output_file.exists():
        print(f"  [OK] Generated structure: {output_file.name}", flush=True)
        return output_file
    else:
        print(f"  [ERROR] Failed to generate structure for {name}", flush=True)
        return None

def generate_orbital_visualization(xyz_file, name):
    """Generate HOMO/LUMO visualizations using PySCF cube files"""
    # Parse XYZ for PySCF
    with open(xyz_file) as f:
        lines = f.readlines()[2:]
    
    atoms_list = []
    for line in lines:
        parts = line.split()
        if len(parts) >= 4:
            atoms_list.append(f"{parts[0]} {float(parts[1]):.6f} {float(parts[2]):.6f} {float(parts[3]):.6f}")
    pyscf_atoms = "; ".join(atoms_list)
    
    homo_file = OUTPUT_DIR / f"mol_{name}_homo.png"
    lumo_file = OUTPUT_DIR / f"mol_{name}_lumo.png"
    
    # Create PySCF script for cube generation
    script = f'''
import sys
import numpy as np
from pyscf import gto, dft
from pyscf.tools import cubegen

try:
    mol = gto.M(
        atom="{pyscf_atoms}",
        basis="3-21g",
        charge=0,
        spin=0,
        symmetry=False,
        verbose=0
    )
    
    mf = dft.RKS(mol, xc="b3lyp")
    mf.kernel()
    
    # Generate cube files
    homo_idx = int(mf.mo_occ.sum() // 2 - 1)  # HOMO index
    lumo_idx = homo_idx + 1  # LUMO index
    
    cubegen.orbital(mol, "{TEMP_DIR}/{name}_homo.cube", mf.mo_coeff[:, homo_idx])
    cubegen.orbital(mol, "{TEMP_DIR}/{name}_lumo.cube", mf.mo_coeff[:, lumo_idx])
    
    print("SUCCESS")
    
except Exception as e:
    print(f"ERROR: {{e}}", file=sys.stderr)
    sys.exit(1)
'''
    
    result = subprocess.run([PYSCF_PYTHON, "-c", script], capture_output=True, text=True, timeout=300)
    
    if "SUCCESS" not in result.stdout:
        print(f"  [ERROR] Failed to generate cubes for {name}: {result.stderr[:100]}", flush=True)
        return None, None
    
    # Fix cube files for xyzrender
    for mo_idx, (cube_file, png_file, orbital) in enumerate([
        (TEMP_DIR / f"{name}_homo.cube", homo_file, "HOMO"),
        (TEMP_DIR / f"{name}_lumo.cube", lumo_file, "LUMO")
    ]):
        if not cube_file.exists():
            continue
        
        # Fix cube file
        lines = open(cube_file).readlines()
        parts = lines[2].split()
        parts[0] = str(-int(parts[0]))
        lines[2] = '  '.join(parts) + '\n'
        natoms = abs(int(parts[0]))
        insert_at = 6 + natoms
        lines.insert(insert_at, f'   1  {mo_idx + 1}\n')
        
        fixed_cube = TEMP_DIR / f"{name}_{orbital.lower()}_fixed.cube"
        with open(fixed_cube, 'w') as f:
            f.writelines(lines)
        
        # Generate visualization
        cmd = f"{XYZRENDER_BIN} {fixed_cube} --mo --flat-mo --iso 0.04 -o {png_file} -t -S 700"
        subprocess.run(cmd, shell=True, capture_output=True, text=True)
        
        if png_file.exists():
            print(f"  [OK] Generated {orbital}: {png_file.name}", flush=True)

def generate_spectra_plot(name, e_s1_ev, osc, e_t1_ev):
    """Generate absorption/emission spectra plot"""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    
    output_file = OUTPUT_DIR / f"mol_{name}_spectra.png"
    
    # Calculate emission
    emission_ev = e_s1_ev - 0.3  # Stokes shift
    emission_nm = 1240.0 / emission_ev
    absorption_nm = 1240.0 / e_s1_ev
    
    # Create plot
    fig, ax = plt.subplots(figsize=(8, 5))
    
    # Simulated absorption spectrum (Gaussian)
    x = np.linspace(200, 800, 500)
    absorption = np.exp(-0.5 * ((1240.0/x - e_s1_ev) / 0.2)**2) * osc
    
    # Simulated emission spectrum
    emission = np.exp(-0.5 * ((1240.0/x - emission_ev) / 0.15)**2) * osc * 0.8
    
    ax.plot(x, absorption, 'b-', label=f'Absorption (λ={absorption_nm:.0f}nm)', linewidth=2)
    ax.plot(x, emission, 'r--', label=f'Emission (λ={emission_nm:.0f}nm)', linewidth=2)
    ax.fill_between(x, absorption, alpha=0.3, color='blue')
    ax.fill_between(x, emission, alpha=0.3, color='red')
    
    ax.set_xlabel('Wavelength (nm)', fontsize=12)
    ax.set_ylabel('Intensity (a.u.)', fontsize=12)
    ax.set_title(f'{name} - Optical Spectra\nΔEST = {e_s1_ev - e_t1_ev:.3f} eV', fontsize=14)
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(300, 700)
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    plt.close()
    
    print(f"  [OK] Generated spectra: {output_file.name}", flush=True)
    return output_file

def main():
    # Load molecules
    with open(DATA_FILE) as f:
        molecules = json.load(f)
    
    print(f"\nLoaded {len(molecules)} molecules from {DATA_FILE.name}", flush=True)
    
    results = []
    skipped_large = []
    
    for name, smiles in molecules.items():
        print(f"\n[{name}]", flush=True)
        print(f"  SMILES: {smiles[:50]}...", flush=True)
        
        # Count heavy atoms
        try:
            heavy_atoms = count_heavy_atoms(smiles)
            print(f"  Heavy atoms: {heavy_atoms}", flush=True)
        except Exception as e:
            print(f"  [ERROR] Failed to count atoms: {e}", flush=True)
            continue
        
        # Skip large molecules
        if heavy_atoms > 35:
            print(f"  [SKIP] Too large for PySCF (>35 heavy atoms)", flush=True)
            skipped_large.append(name)
            continue
        
        # Step 1: Generate 3D conformer
        print("  Generating 3D conformer...", flush=True)
        try:
            xyz_file = generate_3d_conformer(smiles, name)
            if xyz_file is None:
                continue
        except Exception as e:
            print(f"  [ERROR] 3D generation failed: {e}", flush=True)
            continue
        
        # Step 2: xTB optimization
        print("  Running xTB optimization...", flush=True)
        try:
            opt_xyz = run_xtb_optimization(xyz_file, name)
        except Exception as e:
            print(f"  [ERROR] xTB failed: {e}", flush=True)
            opt_xyz = xyz_file
        
        # Step 3: PySCF TDDFT
        print("  Running PySCF TDDFT...", flush=True)
        try:
            e_s1, e_t1, osc = run_pyscf_tddft(opt_xyz, name)
        except Exception as e:
            print(f"  [ERROR] PySCF failed: {e}", flush=True)
            continue
        
        if e_s1 is None:
            continue
        
        # Calculate emission wavelength
        emission_ev = e_s1 - 0.3  # Stokes shift estimate
        emission_nm = 1240.0 / emission_ev
        delta_est = e_s1 - e_t1
        
        # Filter for TADF
        passed = (450 <= emission_nm <= 550) and (delta_est < 0.3) and (osc > 0.001)
        
        result = {
            "name": name,
            "E_S1_ev": round(e_s1, 3),
            "E_T1_ev": round(e_t1, 3),
            "delta_E_ST_ev": round(delta_est, 3),
            "osc": round(osc, 4),
            "emission_nm": round(emission_nm, 1),
            "passed_filter": passed
        }
        results.append(result)
        
        status = "✓ PASSED" if passed else "✗ FAILED"
        print(f"  E_S1={e_s1:.3f}eV, E_T1={e_t1:.3f}eV, ΔEST={delta_est:.3f}eV", flush=True)
        print(f"  Emission={emission_nm:.1f}nm, Osc={osc:.4f} → {status}", flush=True)
    
    # Step 4: Visualize top candidates
    print("\n" + "=" * 60, flush=True)
    print("Generating Visualizations for Top Candidates", flush=True)
    print("=" * 60, flush=True)
    
    if not results:
        print("\n[ERROR] No molecules processed successfully!", flush=True)
        return
    
    # Sort by ΔEST (ascending) and filter
    top_candidates = sorted([r for r in results if r["passed_filter"]], 
                           key=lambda x: x["delta_E_ST_ev"])[:5]
    
    if not top_candidates:
        print("No molecules passed the filter. Visualizing top 3 by ΔEST...", flush=True)
        top_candidates = sorted(results, key=lambda x: x["delta_E_ST_ev"])[:3]
    
    print(f"\nTop {len(top_candidates)} candidates:", flush=True)
    for i, cand in enumerate(top_candidates, 1):
        print(f"  {i}. {cand['name']}: ΔEST={cand['delta_E_ST_ev']:.3f}eV, λ={cand['emission_nm']:.0f}nm", flush=True)
    
    for candidate in top_candidates:
        name = candidate["name"]
        print(f"\n[{name}] Generating visualizations...", flush=True)
        
        xyz_file = TEMP_DIR / f"{name}_xtbopt.xyz"
        if not xyz_file.exists():
            xyz_file = TEMP_DIR / f"{name}_initial.xyz"
        
        if not xyz_file.exists():
            print(f"  [ERROR] XYZ file not found", flush=True)
            continue
        
        # Structure
        try:
            generate_structure_visualization(xyz_file, name)
        except Exception as e:
            print(f"  [ERROR] Structure viz failed: {e}", flush=True)
        
        # HOMO/LUMO
        try:
            generate_orbital_visualization(xyz_file, name)
        except Exception as e:
            print(f"  [ERROR] Orbital viz failed: {e}", flush=True)
        
        # Spectra
        try:
            generate_spectra_plot(name, candidate["E_S1_ev"], candidate["osc"], candidate["E_T1_ev"])
        except Exception as e:
            print(f"  [ERROR] Spectra viz failed: {e}", flush=True)
    
    # Save results to CSV
    csv_file = OUTPUT_DIR / "screening_results.csv"
    with open(csv_file, 'w') as f:
        f.write("name,E_S1_ev,E_T1_ev,delta_E_ST_ev,osc,emission_nm,passed_filter\n")
        for r in results:
            f.write(f"{r['name']},{r['E_S1_ev']},{r['E_T1_ev']},{r['delta_E_ST_ev']},"
                   f"{r['osc']},{r['emission_nm']},{r['passed_filter']}\n")
    
    print(f"\n" + "=" * 60, flush=True)
    print(f"Results saved to: {csv_file}", flush=True)
    print(f"Visualizations saved to: {OUTPUT_DIR}", flush=True)
    print("=" * 60, flush=True)
    
    # Summary
    passed_count = sum(1 for r in results if r["passed_filter"])
    print(f"\nSummary:", flush=True)
    print(f"  Processed: {len(results)}/{len(molecules)} molecules", flush=True)
    print(f"  Passed TADF filter: {passed_count}/{len(results)}", flush=True)
    print(f"  Skipped (too large): {len(skipped_large)} - {', '.join(skipped_large)}", flush=True)
    print(f"  Filter criteria: 450-550nm emission, ΔEST < 0.3eV, osc > 0.001", flush=True)

if __name__ == "__main__":
    main()
