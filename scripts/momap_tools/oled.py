#!/usr/bin/env python3
"""
momap-oled — Extended OLED tool wrappers for isc_tvcf, ic_tvcf, pysoc, sumstat, transport.
"""
import os
import sys
import subprocess
import argparse
from pathlib import Path

TOOLS_DIR = os.path.dirname(os.path.abspath(__file__))

# ─── isc_tvcf ───────────────────────────────────────────────────────────
def generate_isc_input(evc_dat, ead_au, hso_cm1, output='momap_isc.inp', temp=298, tmax=5000):
    """Generate isc_tvcf input for phosphorescence spectrum + k_ISC."""
    content = f"""do_isc_tvcf_ft   = 1
do_isc_tvcf_spec = 1

&isc_tvcf
  DUSHIN  = .t.
  Temp    = {temp} K
  tmax    = {tmax} fs
  dt      = 0.001 fs
  Ead     = {ead_au:.8f} au
  Hso     = {hso_cm1:.6f} cm-1
  DSFile  = "{evc_dat}"
  Emax    = 0.3 au
  logFile = "isc.tvcf.log"
  FtFile  = "isc.tvcf.ft.dat"
  FoFile  = "isc.tvcf.fo.dat"
/
"""
    with open(output, 'w') as f:
        f.write(content)
    return output


def parse_isc_output(fo_file='isc.tvcf.fo.dat'):
    """Parse isc_tvcf output for k_ISC and k_RISC."""
    with open(fo_file) as f:
        text = f.read()
    import re
    k_isc = re.search(r'rate is\s+([\d.E+-]+)\s+s-1', text.split('\n')[0])
    k_risc = re.search(r'rate is\s+([\d.E+-]+)\s+s-1', text.split('\n')[1])
    return {
        'k_ISC_s1': float(k_isc.group(1)) if k_isc else None,
        'k_RISC_s1': float(k_risc.group(1)) if k_risc else None,
    }


# ─── ic_tvcf ────────────────────────────────────────────────────────────
def generate_ic_input(evc_dat, nac_coul, ead_au, output='momap_ic.inp', temp=300, tmax=655):
    """Generate ic_tvcf input for internal conversion rate."""
    content = f"""do_ic_tvcf_ft   = 1
do_ic_tvcf_spec = 1

&ic_tvcf
  DUSHIN   = .t.
  Temp     = {temp} K
  tmax     = {tmax} fs
  dt       = 0.01 fs
  Ead      = {ead_au:.8f} au
  DSFile   = "{evc_dat}"
  CoulFile = "{nac_coul}"
  logFile  = "ic.tvcf.log"
  FtFile   = "ic.tvcf.ft.dat"
  FoFile   = "ic.tvcf.fo.dat"
/
"""
    with open(output, 'w') as f:
        f.write(content)
    return output


# ─── pysoc ──────────────────────────────────────────────────────────────
def generate_pysoc_input(qm_input_com, output='momap_pysoc.inp',
                         n_singlets=4, n_triplets=4,
                         qc_exe='g16', qc_ppn=8):
    """Generate PySOC input for spin-orbit coupling calculation."""
    content = f"""do_pysoc = 1

&pysoc
  sched_type = local
  qc_exe     = {qc_exe}
  qc_ppn     = {qc_ppn}
  pysoc_QM_code = 'gauss_tddft'
  pysoc_QM_input_file = {qm_input_com}
  n_excited_singlets = {n_singlets}
  n_excited_triplets = {n_triplets}
/
"""
    with open(output, 'w') as f:
        f.write(content)
    return output


# ─── spec_sums ──────────────────────────────────────────────────────────
def generate_sums_input(evc_dat, ead_au, dipole_abs, dipole_emi,
                        output='momap_sums.inp', fwhm=500, maxvib=10):
    """Generate spec_sums input — all rates + quantum yield in one pass."""
    content = f"""do_spec_sums = 1

&spec_sums
  DSFile     = "{evc_dat}"
  Ead        = {ead_au:.8f} au
  dipole_abs = {dipole_abs:.6f} debye
  dipole_emi = {dipole_emi:.6f} debye
  maxvib     = {maxvib}
  if_cal_ic  = .t.
  FWHM       = {fwhm} cm-1
  flog       = "spec.sums.log"
/
"""
    with open(output, 'w') as f:
        f.write(content)
    return output


def parse_sums_output(logfile='spec.sums.log'):
    """Parse spec_sums output for all rates."""
    with open(logfile) as f:
        text = f.read()
    import re
    k_emi = re.search(r'Emission rate is\s+([\d.E+-]+)\s+s-1', text)
    return {
        'k_r_s1': float(k_emi.group(1)) if k_emi else None,
    }


# ─── transport ──────────────────────────────────────────────────────────
def generate_transport_input(cif_file, output='momap_transport.inp',
                             basis='b3lyp STO-3g', temp=300,
                             qc_exe='g16', ratetype='marcus',
                             nsimu=2000, tsimu=1000):
    """Generate transport input for charge carrier mobility."""
    content = f"""&transport
  do_transport_prepare              = 1
  do_transport_submit_HL_job        = 1
  do_transport_get_transferintegral = 1
  do_transport_submit_RE_job        = 1
  do_transport_get_re_evc           = 1
  do_transport_run_MC               = 1
  do_transport_get_mob_MC           = 1
  do_transport_gather_momap_data    = 1

  sched_type     = local
  queue_name     = localhost
  compute_engine = 1
  qc_exe         = {qc_exe}
  basis_name     = {basis}
  basis_name_re  = {basis}
  qc_memory      = 4096
  qc_nodes       = 1
  qc_ppn         = 2
  temp           = {temp}
  ratetype       = {ratetype}
  lat_cutoff     = 4
  nsimu          = {nsimu}
  tsimu          = {tsimu}
  crystal        = {cif_file}
/
"""
    with open(output, 'w') as f:
        f.write(content)
    return output


# ─── CLI ────────────────────────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser(description='MOMAP OLED tool wrappers')
    sub = parser.add_subparsers(dest='cmd')

    p = sub.add_parser('isc', help='Generate isc_tvcf input')
    p.add_argument('--evc-dat', default='evc.cart.dat')
    p.add_argument('--ead', type=float, required=True)
    p.add_argument('--hso', type=float, default=1.0)
    p.add_argument('-o', default='momap_isc.inp')

    p = sub.add_parser('ic', help='Generate ic_tvcf input')
    p.add_argument('--evc-dat', default='evc.cart.dat')
    p.add_argument('--nac', default='evc.cart.nac')
    p.add_argument('--ead', type=float, required=True)
    p.add_argument('-o', default='momap_ic.inp')

    p = sub.add_parser('pysoc', help='Generate PySOC input')
    p.add_argument('--com', required=True, help='Gaussian TDDFT input file')
    p.add_argument('--n-singlets', type=int, default=4)
    p.add_argument('--n-triplets', type=int, default=4)
    p.add_argument('-o', default='momap_pysoc.inp')

    p = sub.add_parser('sums', help='Generate spec_sums input')
    p.add_argument('--evc-dat', default='evc.cart.dat')
    p.add_argument('--ead', type=float, required=True)
    p.add_argument('--dipole-abs', type=float, required=True)
    p.add_argument('--dipole-emi', type=float, required=True)
    p.add_argument('-o', default='momap_sums.inp')

    p = sub.add_parser('transport', help='Generate transport input')
    p.add_argument('--cif', required=True, help='Crystal structure .cif file')
    p.add_argument('-o', default='momap_transport.inp')

    args = parser.parse_args()
    if not args.cmd:
        parser.print_help()
        return 1

    if args.cmd == 'isc':
        path = generate_isc_input(args.evc_dat, args.ead, args.hso, args.o)
        print(f"✅ isc_tvcf input → {path}")
    elif args.cmd == 'ic':
        path = generate_ic_input(args.evc_dat, args.nac, args.ead, args.o)
        print(f"✅ ic_tvcf input → {path}")
    elif args.cmd == 'pysoc':
        path = generate_pysoc_input(args.com, args.o, args.n_singlets, args.n_triplets)
        print(f"✅ pysoc input → {path}")
    elif args.cmd == 'sums':
        path = generate_sums_input(args.evc_dat, args.ead, args.dipole_abs, args.dipole_emi, args.o)
        print(f"✅ spec_sums input → {path}")
    elif args.cmd == 'transport':
        path = generate_transport_input(args.cif, args.o)
        print(f"✅ transport input → {path}")

    return 0

if __name__ == '__main__':
    sys.exit(main())
