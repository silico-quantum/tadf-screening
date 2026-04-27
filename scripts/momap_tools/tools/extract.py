#!/usr/bin/env python3
"""
momap-extract — Parse Gaussian output files, extract MOMAP parameters.
Generates a complete spec_tvcf input file ready to run.
"""
import sys
import re
import os
import argparse
from pathlib import Path

AU2DEBYE = 2.541746

def extract_scf_energy(logpath):
    """Extract last SCF Done energy from Gaussian log."""
    energies = []
    with open(logpath) as f:
        for line in f:
            m = re.search(r'SCF Done:\s+E\([^)]+\)\s*=\s*(-?\d+\.\d+)', line)
            if m:
                energies.append(float(m.group(1)))
    if not energies:
        raise ValueError(f"No SCF Done found in {logpath}")
    return energies[-1]

def extract_last_excitation_ev(logpath):
    """Extract the last S1 (state 1) excitation energy (eV) from TDDFT log.
    Returns the adiabatic excitation energy = the very last Excited State 1 value."""
    last_exc_ev = None
    with open(logpath) as f:
        for line in f:
            # Match: " Excited State   1:      Singlet-A'     1.6143 eV  768.02 nm  f=0.0026"
            if 'Excited State' in line and ':      ' in line:
                parts = line.split()
                try:
                    state_idx = parts.index('State') if 'State' in parts else -1
                    if state_idx >= 0:
                        state_num = int(parts[state_idx + 1].rstrip(':'))
                        if state_num == 1:
                            # Find eV value
                            for j, p in enumerate(parts):
                                if p == 'eV':
                                    last_exc_ev = float(parts[j-1])
                                    break
                except (ValueError, IndexError):
                    pass
    return last_exc_ev

def extract_transition_dipoles(logpath):
    """Extract transition electric dipole moments for all excited states.
    Returns list of dicts with state, X, Y, Z, Dip.S., Osc."""
    results = []
    with open(logpath) as f:
        lines = f.readlines()
    
    in_block = False
    for line in lines:
        if 'Ground to excited state transition electric dipole moments' in line:
            in_block = True
            continue
        if in_block:
            m = re.match(r'\s+(\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)', line)
            if m:
                results.append({
                    'state': int(m.group(1)),
                    'X': float(m.group(2)),
                    'Y': float(m.group(3)),
                    'Z': float(m.group(4)),
                    'DipS': float(m.group(5)),
                    'Osc': float(m.group(6)),
                })
            elif line.strip() == '' or 'Ground to excited state transition velocity' in line:
                in_block = False
    return results

def extract_dipole_moment(logpath):
    """Extract last dipole moment from Gaussian log."""
    dipoles = []
    with open(logpath) as f:
        for line in f:
            m = re.match(r'\s+X=\s*(-?\d+\.\d+)\s+Y=\s*(-?\d+\.\d+)\s+Z=\s*(-?\d+\.\d+)\s+Tot=\s*(-?\d+\.\d+)', line)
            if m:
                dipoles.append({
                    'X': float(m.group(1)), 'Y': float(m.group(2)),
                    'Z': float(m.group(3)), 'Tot': float(m.group(4))
                })
    return dipoles if dipoles else None

def count_normal_terminations(logpath):
    """Count Normal termination occurrences."""
    count = 0
    with open(logpath) as f:
        for line in f:
            if 'Normal termination' in line:
                count += 1
    return count

def generate_spec_tvcf_input(params, output_path='momap_spec.inp'):
    """Generate MOMAP spec_tvcf input file."""
    template = f"""do_spec_tvcf_ft   = 1
do_spec_tvcf_spec = 1

&spec_tvcf
  DUSHIN   = .t.
  HERZ     = .t.
  Temp     = {params.get('temperature', 300)} K
  tmax     = {params.get('tmax', 5000)} fs
  dt       = {params.get('dt', 0.001)} fs
  Ead      = {params.get('Ead', 0.07509):.8f} au
  EDMA     = {params.get('EDMA', 0.92694):.6f} debye
  EDME     = {params.get('EDME', 0.64751):.6f} debye
  FreqScale  = 1.0
  DSFile     = "{params.get('dsfile', 'evc.cart.dat')}"
  Emax       = {params.get('emax', 0.3)} au
  dE         = {params.get('de', 0.00001)} au
  logFile    = "spec.tvcf.log"
  FtFile     = "spec.tvcf.ft.dat"
  FoFile     = "spec.tvcf.fo.dat"
  FoSFile    = "spec.tvcf.spec.dat"
/
"""
    with open(output_path, 'w') as f:
        f.write(template)
    return output_path

def main():
    parser = argparse.ArgumentParser(description='Extract MOMAP parameters from Gaussian logs')
    parser.add_argument('--s0', required=True, help='S0 ground state Gaussian log')
    parser.add_argument('--s1', required=True, help='S1 excited state Gaussian log')
    parser.add_argument('--t1', help='T1 triplet state Gaussian log (optional)')
    parser.add_argument('--output', '-o', default='momap_spec.inp', help='Output MOMAP input file')
    parser.add_argument('--temperature', '-T', type=float, default=300, help='Temperature (K)')
    parser.add_argument('--json', action='store_true', help='Output JSON instead of input file')
    import json as json_mod

    args = parser.parse_args()

    # Extract SCF energies + adiabatic excitation
    E_S0 = extract_scf_energy(args.s0)
    E_S1_scf = extract_scf_energy(args.s1)  # S0 part at S1 geometry
    exc_ev = extract_last_excitation_ev(args.s1)
    
    if exc_ev:
        # Ead = E(S1 at S1_min) - E(S0 at S0_min)
        #      = (E_SCF_S1 + E_exc) - E_SCF_S0
        E_S1_total = E_S1_scf + exc_ev / 27.2114
        Ead = E_S1_total - E_S0
    else:
        Ead = E_S1_scf - E_S0  # fallback

    # Extract transition dipoles from S1 log
    # EDMA: from first TDM block (S0 geometry → vertical absorption)
    # EDME: from last TDM block (S1 minimum → emission)
    tdms = extract_transition_dipoles(args.s1)
    
    EDMA_au = tdms[0]['DipS'] if tdms else 0.0
    EDME_au = tdms[-1]['DipS'] if tdms else 0.0
    EDMA = EDMA_au * AU2DEBYE
    EDME = EDME_au * AU2DEBYE

    # Extract oscillator strengths
    f_abs = tdms[0]['Osc'] if tdms else 0.0
    f_emi = tdms[-1]['Osc'] if tdms else 0.0

    # Count normal terminations
    nts_s0 = count_normal_terminations(args.s0)
    nts_s1 = count_normal_terminations(args.s1)
    nts_t1 = count_normal_terminations(args.t1) if args.t1 else None

    params = {
        'Ead': Ead,
        'EDMA': EDMA,
        'EDME': EDME,
        'f_abs': f_abs,
        'f_emi': f_emi,
        'E_S0': E_S0,
        'E_S1': E_S1_total if exc_ev else E_S1_scf,
        'E_S1_scf': E_S1_scf,
        'E_exc_ev': exc_ev,
        'temperature': args.temperature,
        'dsfile': 'evc.cart.dat',
        'tmax': 5000,
        'dt': 0.001,
        'emax': 0.3,
        'de': 0.00001,
        'nts_s0': nts_s0,
        'nts_s1': nts_s1,
        'nts_t1': nts_t1,
    }

    if args.json:
        print(json_mod.dumps(params, indent=2))
    else:
        path = generate_spec_tvcf_input(params, args.output)
        print(f"✅ MOMAP spec_tvcf input → {path}")
        print(f"   Ead  = {Ead:.8f} au ({Ead*27.2114:.4f} eV)")
        print(f"   EDMA = {EDMA:.4f} debye (f={f_abs:.4f})")
        print(f"   EDME = {EDME:.4f} debye (f={f_emi:.4f})")
        print(f"   S0 NTs: {nts_s0}, S1 NTs: {nts_s1}" + (f", T1 NTs: {nts_t1}" if nts_t1 else ""))

if __name__ == '__main__':
    main()
