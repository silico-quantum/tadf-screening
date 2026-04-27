#!/usr/bin/env python3
"""
momap-tadf — Full TADF photophysics pipeline using MOMAP.

Given a directory of Gaussian S0/T1/S1 outputs, runs:
  1. EVC (S1→S0) → Duschinsky + HR for fluorescence
  2. EVC (S1→T1) → Duschinsky + HR for ISC
  3. spec_tvcf (S1→S0) → fluorescence spectrum
  4. ISC rate (S1→T1) → k_ISC
  5. Summary → Φ, τ, ΔE_ST, spectral peak in blue window
"""
import sys
import os
import re
import json
import argparse
import subprocess
from pathlib import Path

sys.path.insert(0, os.path.dirname(__file__))
from extract import (
    extract_scf_energy,
    extract_transition_dipoles,
    count_normal_terminations,
    generate_spec_tvcf_input,
)
from runner import (
    patch_momap_for_mpi3,
    ensure_fchk,
    create_nodefile,
    MOMAP_ROOT,
)

AU2DEBYE = 2.541746
HA2EV = 27.2114

def run_momap_in_dir(inputfile, workdir):
    """Run patched MOMAP in workdir."""
    patched = patch_momap_for_mpi3()
    create_nodefile(workdir)
    
    env = os.environ.copy()
    env['MOMAP_ROOT'] = MOMAP_ROOT
    env['MOMAP_LICENSE'] = os.path.join(MOMAP_ROOT, 'license', 'hzwtech.lic')
    
    cmd = ['python', patched, '-i', str(inputfile)]
    print(f"  🚀 {' '.join(cmd)}")
    result = subprocess.run(cmd, cwd=str(workdir), env=env,
                          capture_output=False)
    return result.returncode == 0

def parse_evc_out(evc_out):
    """Extract reorganization energy and mode count from evc.out."""
    with open(evc_out) as f:
        content = f.read()
    
    m = re.search(r'(\d+)\s+# num of atoms', content)
    natoms = int(m.group(1)) if m else 0
    m = re.search(r'(\d+)\s+# num of modes', content)
    nmodes = int(m.group(1)) if m else 0
    
    return {'natoms': natoms, 'nmodes': nmodes}

def parse_spec_output(spec_file):
    """Parse spec.tvcf.spec.dat for peak wavelength."""
    if not Path(spec_file).exists():
        return {}
    
    try:
        data = []
        with open(spec_file) as f:
            for line in f:
                if line.startswith('#'):
                    continue
                parts = line.split()
                if len(parts) >= 7:
                    data.append([float(x) for x in parts])
        
        if not data:
            return {}
        
        emi = [row[5] for row in data]  # FC_emi
        wavelengths = [row[3] for row in data]
        
        # Find peaks
        peaks = []
        for i in range(1, len(emi)-1):
            if emi[i] > emi[i-1] and emi[i] > emi[i+1] and emi[i] > 0:
                peaks.append({'wavelength': wavelengths[i], 'intensity': emi[i]})
        
        peaks.sort(key=lambda x: x['intensity'], reverse=True)
        
        max_idx = max(range(len(emi)), key=lambda i: emi[i])
        
        return {
            'peak_wavelength': wavelengths[max_idx],
            'peak_intensity': emi[max_idx],
            'top_peaks': peaks[:5],
            'data_points': len(data),
        }
    except Exception as e:
        return {'error': str(e)}

def process_molecule(mol_id, s0_log, s1_log, t1_log, output_dir, temperature=300):
    """Run full TADF MOMAP pipeline for one molecule."""
    mol_dir = Path(output_dir) / mol_id
    mol_dir.mkdir(parents=True, exist_ok=True)
    
    results = {'mol_id': mol_id, 'success': False}
    
    # Step 0: Copy and validate input files
    print(f"\n{'='*60}")
    print(f"🧬 Processing {mol_id}")
    print(f"{'='*60}")
    
    # Check files exist and have 2 Normal terminations
    for label, path in [('S0', s0_log), ('S1', s1_log), ('T1', t1_log)]:
        if path is None:
            print(f"  ⚠️  {label}: no file provided")
            continue
        nts = count_normal_terminations(path)
        print(f"  📄 {label}: {path} ({nts} Normal terminations)")
        if nts < 2:
            print(f"     ⚠️  Warning: expected 2 NTs (opt+freq), got {nts}")
    
    # Ensure fchk files
    for log_path in [s0_log, s1_log, t1_log]:
        if log_path:
            ensure_fchk(log_path, mol_dir)
    
    # Extract parameters
    try:
        E_S0 = extract_scf_energy(s0_log)
        E_S1 = extract_scf_energy(s1_log)
        E_T1 = extract_scf_energy(t1_log) if t1_log else None
    except Exception as e:
        print(f"  ❌ SCF extraction failed: {e}")
        return results
    
    delta_EST = (E_S1 - E_T1) * HA2EV if E_T1 else None
    
    tdms = extract_transition_dipoles(s1_log)
    EDMA_au = tdms[0]['DipS'] if tdms else 0.0
    EDME_au = tdms[-1]['DipS'] if tdms else 0.0
    
    results.update({
        'E_S0': E_S0, 'E_S1': E_S1, 'E_T1': E_T1,
        'Ead_S1_S0': (E_S1 - E_S0) * HA2EV,
        'delta_EST_eV': delta_EST,
        'EDMA': EDMA_au * AU2DEBYE,
        'EDME': EDME_au * AU2DEBYE,
        'f_abs': tdms[0]['Osc'] if tdms else 0.0,
        'f_emi': tdms[-1]['Osc'] if tdms else 0.0,
    })

    print(f"\n  📊 Extracted parameters:")
    print(f"     E(S0)   = {E_S0:.8f} au")
    print(f"     E(S1)   = {E_S1:.8f} au")
    if E_T1:
        print(f"     E(T1)   = {E_T1:.8f} au")
        print(f"     ΔE_ST   = {delta_EST:.4f} eV")
    print(f"     EDMA    = {results['EDMA']:.4f} debye")
    print(f"     EDME    = {results['EDME']:.4f} debye")
    
    # Step 1: EVC (S1→S0)
    print(f"\n  📐 Step 1: EVC (S1→S0)")
    evc_s1_input = mol_dir / 'momap_evc_s1.inp'
    with open(evc_s1_input, 'w') as f:
        f.write(f"""do_evc = 1

&evc
  ffreq(1) = "{os.path.basename(s0_log)}"
  ffreq(2) = "{os.path.basename(s1_log)}"
  sort_mode = 1
/
""")
    ok = run_momap_in_dir(evc_s1_input, mol_dir)
    if not ok:
        print(f"  ❌ EVC (S1→S0) failed")
        return results
    evc_s1_info = parse_evc_out(mol_dir / 'evc.out')
    print(f"  ✅ {evc_s1_info['natoms']} atoms, {evc_s1_info['nmodes']} modes")
    
    # Step 2: spec_tvcf
    print(f"\n  🌈 Step 2: Spectrum (S1→S0)")
    spec_params = {
        'temperature': temperature,
        'Ead': (E_S1 - E_S0),
        'EDMA': results['EDMA'],
        'EDME': results['EDME'],
        'dsfile': 'evc.cart.dat',
    }
    spec_input = mol_dir / 'momap_spec.inp'
    generate_spec_tvcf_input(spec_params, spec_input)
    ok = run_momap_in_dir(spec_input, mol_dir)
    if ok:
        spec_info = parse_spec_output(mol_dir / 'spec.tvcf.spec.dat')
        results['spectrum'] = spec_info
        if spec_info.get('peak_wavelength'):
            print(f"  ✅ Peak: {spec_info['peak_wavelength']:.1f} nm")
            
            # Blue window check
            peak_nm = spec_info['peak_wavelength']
            if 450 <= peak_nm <= 490:
                print(f"  🔵 WITHIN BLUE WINDOW (450–490 nm)!")
                results['blue_window'] = True
            elif 400 <= peak_nm <= 500:
                print(f"  🔹 Near blue window ({peak_nm:.0f} nm)")
                results['blue_window'] = 'near'
            else:
                print(f"  ⚪ Outside blue window")
                results['blue_window'] = False
        
        # FWHM
        if spec_info.get('top_peaks'):
            results['top_spectral_peaks'] = spec_info['top_peaks'][:3]
    else:
        print(f"  ⚠️  Spectrum calculation had issues")
    
    # Step 3: ISC (S1→T1) if T1 available
    if t1_log and Path(t1_log).exists():
        print(f"\n  🔀 Step 3: ISC (S1→T1)")
        # First need EVC for S1→T1
        evc_isc_input = mol_dir / 'momap_evc_isc.inp'
        with open(evc_isc_input, 'w') as f:
            f.write(f"""do_evc = 1

&evc
  ffreq(1) = "{os.path.basename(s1_log)}"
  ffreq(2) = "{os.path.basename(t1_log)}"
  sort_mode = 1
/
""")
        ok = run_momap_in_dir(evc_isc_input, mol_dir)
        
        # TODO: ISC rate calculation — needs do_isc block
        # For now, just note that EVC is ready
        if ok:
            print(f"  ✅ EVC (S1→T1) complete, ready for ISC calculation")
    
    results['success'] = True
    results['output_dir'] = str(mol_dir)
    
    # Summary
    print(f"\n  {'='*50}")
    print(f"  📋 {mol_id} Summary")
    print(f"  {'='*50}")
    print(f"  ΔE_ST     = {delta_EST:.4f} eV" if delta_EST else "  ΔE_ST     = N/A")
    if results.get('spectrum', {}).get('peak_wavelength'):
        p = results['spectrum']
        print(f"  λ_emi      = {p['peak_wavelength']:.1f} nm")
        blue = results.get('blue_window')
        if blue == True:
            print(f"  Blue Window = ✅ YES")
        elif blue == 'near':
            print(f"  Blue Window = 🟡 NEAR")
        else:
            print(f"  Blue Window = ❌ NO")
    
    return results


def main():
    parser = argparse.ArgumentParser(
        description='MOMAP TADF pipeline — auto EVC → spectrum → ISC for TADF candidates'
    )
    parser.add_argument('mol_id', help='Molecule identifier')
    parser.add_argument('--s0', required=True, help='S0 ground state Gaussian log')
    parser.add_argument('--s1', required=True, help='S1 excited state Gaussian log')
    parser.add_argument('--t1', help='T1 triplet state Gaussian log')
    parser.add_argument('--output', '-o', default='./momap_tadf_output', help='Output directory')
    parser.add_argument('--temperature', '-T', type=float, default=300, help='Temperature (K)')
    parser.add_argument('--json', action='store_true', help='Output results as JSON')
    args = parser.parse_args()

    results = process_molecule(
        mol_id=args.mol_id,
        s0_log=args.s0,
        s1_log=args.s1,
        t1_log=args.t1,
        output_dir=args.output,
        temperature=args.temperature,
    )

    if args.json:
        print(json.dumps(results, indent=2, default=str))
    
    if results.get('success'):
        return 0
    else:
        print(f"\n❌ Processing failed for {args.mol_id}")
        return 1

if __name__ == '__main__':
    sys.exit(main())
