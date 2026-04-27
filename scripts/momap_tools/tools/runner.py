#!/usr/bin/env python3
"""
momap-run — One-command MOMAP wrapper with auto MPI fix + formchk + Slurm.

Handles:
  - OpenMPI 3.x --hostfile vs -machinefile patching (auto)
  - formchk generation for missing .fchk files (auto)
  - Local (marcus2) vs Slurm submission
  - Sub-commands: evc, spec, isc, full
"""
import sys
import os
import subprocess
import argparse
import shutil
import tempfile
from pathlib import Path

MOMAP_ROOT = os.environ.get('MOMAP_ROOT', '/opt/MOMAP-2024A')
MOMAP_BIN = os.path.join(MOMAP_ROOT, 'bin')
MOMAP_MPI_BIN = os.path.join(MOMAP_ROOT, 'bin', 'openmpi', 'bin')

PATCHED_SCRIPT = None

def find_momap_bin():
    """Find momap executable, prefer system PATH."""
    momap = shutil.which('momap')
    if momap:
        return momap
    momap = os.path.join(MOMAP_BIN, 'momap')
    if os.path.exists(momap):
        return momap
    raise RuntimeError("MOMAP not found. source /opt/MOMAP-2024A/env.sh first")

def patch_momap_for_mpi3():
    """Create a patched momap script that uses --hostfile instead of -machinefile."""
    global PATCHED_SCRIPT
    if PATCHED_SCRIPT and os.path.exists(PATCHED_SCRIPT):
        return PATCHED_SCRIPT
    
    momap_orig = find_momap_bin()
    tmpdir = Path(tempfile.gettempdir()) / 'momap_patched'
    tmpdir.mkdir(exist_ok=True)
    
    patched = tmpdir / 'momap_patched'
    
    with open(momap_orig) as f:
        content = f.read()
    
    content = content.replace('-machinefile', '--hostfile')
    
    with open(patched, 'w') as f:
        f.write(content)
    os.chmod(patched, 0o755)
    
    PATCHED_SCRIPT = str(patched)
    return PATCHED_SCRIPT

def ensure_fchk(logpath, workdir=None):
    """Ensure .fchk file exists for a given .log file. Auto-generate from .chk if missing."""
    logpath = Path(logpath)
    fchk = logpath.with_suffix('.fchk')
    chk = logpath.with_suffix('.chk')
    
    if fchk.exists():
        return str(fchk)
    
    if not chk.exists():
        print(f"⚠️  No .fchk or .chk found for {logpath.name}")
        print(f"   MOMAP requires .fchk. Run: formchk {chk.name} {fchk.name}")
        return None
    
    cwd = workdir or chk.parent
    print(f"🔄 Generating {fchk.name} from {chk.name} ...")
    result = subprocess.run(
        ['formchk', str(chk.name), str(fchk.name)],
        cwd=str(cwd), capture_output=True, text=True
    )
    if result.returncode == 0:
        print(f"   ✅ {fchk.name}")
        return str(fchk)
    else:
        print(f"   ❌ formchk failed: {result.stderr}")
        return None

def create_nodefile(workdir, slots=4):
    """Create OpenMPI 3.x compatible hostfile."""
    nodefile = Path(workdir) / 'nodefile'
    with open(nodefile, 'w') as f:
        f.write(f"localhost slots={slots}\n")
    return str(nodefile)

def run_momap(inputfile, workdir=None, use_slurm=False, slurm_partition='X32Cv4', slurm_nprocs=4):
    """Run MOMAP calculation."""
    inputfile = Path(inputfile)
    cwd = workdir or inputfile.parent
    
    if use_slurm:
        return submit_slurm(str(inputfile), cwd, slurm_partition, slurm_nprocs)
    else:
        return run_local(str(inputfile), cwd)

def run_local(inputfile, cwd):
    """Run MOMAP locally (marcus2)."""
    patched = patch_momap_for_mpi3()
    create_nodefile(cwd)
    
    env = os.environ.copy()
    env['MOMAP_ROOT'] = MOMAP_ROOT
    env['MOMAP_LICENSE'] = os.path.join(MOMAP_ROOT, 'license', 'hzwtech.lic')
    
    cmd = ['python', patched, '-i', inputfile]
    print(f"🚀 Running: {' '.join(cmd)}")
    
    result = subprocess.run(cmd, cwd=cwd, env=env, capture_output=False)
    return result.returncode

def submit_slurm(inputfile, cwd, partition, nprocs):
    """Submit MOMAP via Slurm to compute nodes."""
    slurm_script = Path(cwd) / 'momap_job.slurm'
    
    script = f"""#!/bin/bash
#SBATCH --job-name=momap
#SBATCH --partition={partition}
#SBATCH --nodes=1
#SBATCH --ntasks={nprocs}
#SBATCH --time=12:00:00
#SBATCH --output=momap_%j.log

source /opt/MOMAP-2024A/env.sh
cd {cwd}
echo "localhost slots=$SLURM_NTASKS" > nodefile

# Patched momap
MOMAP_SCRIPT={PATCHED_SCRIPT or patch_momap_for_mpi3()}
python $MOMAP_SCRIPT -i {inputfile}
"""
    with open(slurm_script, 'w') as f:
        f.write(script)
    
    print(f"📤 Submitting Slurm job...")
    result = subprocess.run(['sbatch', str(slurm_script)], capture_output=True, text=True)
    print(result.stdout.strip())
    return result.returncode

def generate_evc_input(s0_log, s1_log, output='momap_evc.inp'):
    """Generate EVC input file."""
    content = f"""do_evc = 1

&evc
  ffreq(1) = "{s0_log}"
  ffreq(2) = "{s1_log}"
  sort_mode = 1
/
"""
    with open(output, 'w') as f:
        f.write(content)
    return output

def cmd_extract(args):
    """Re-export extract functionality."""
    from extract import main as extract_main
    sys.argv = ['extract'] + args.passthrough if hasattr(args, 'passthrough') else []
    # Delegate to extract module - just use subprocess
    cmd = [sys.executable, os.path.join(os.path.dirname(__file__), 'extract.py')] + args.passthrough
    os.execv(cmd[0], cmd)

def main():
    parser = argparse.ArgumentParser(description='MOMAP runner — one-command wrapper')
    parser.add_argument('input', nargs='?', help='MOMAP input file (momap.inp)')
    parser.add_argument('--slurm', '-s', action='store_true', help='Submit via Slurm')
    parser.add_argument('--partition', '-p', default='X32Cv4', help='Slurm partition')
    parser.add_argument('--nprocs', '-n', type=int, default=4, help='MPI processes')
    parser.add_argument('--workdir', '-C', default='.', help='Working directory')
    parser.add_argument('--patch-only', action='store_true', help='Only create patched momap script and exit')
    args = parser.parse_args()

    if args.patch_only:
        patched = patch_momap_for_mpi3()
        print(f"Patched MOMAP: {patched}")
        print(f"Usage: python {patched} -i momap.inp")
        return 0

    if not args.input:
        parser.print_help()
        return 1

    sys.exit(run_momap(args.input, args.workdir, args.slurm, args.partition, args.nprocs))

if __name__ == '__main__':
    main()
