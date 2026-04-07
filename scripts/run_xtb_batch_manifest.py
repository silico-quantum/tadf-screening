#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import re
import subprocess
from pathlib import Path

ENERGY_RE = re.compile(r"TOTAL ENERGY\s+(-?\d+\.\d+)\s+Eh")
GAP_RE = re.compile(r"HOMO-LUMO GAP\s+(-?\d+\.\d+)\s+eV")


def parse_xtb_out(out_path: Path) -> dict:
    text = out_path.read_text(errors='ignore') if out_path.exists() else ''
    m_e = ENERGY_RE.search(text)
    m_g = GAP_RE.search(text)
    ok_term = 'normal termination of xtb' in text
    return {
        'normal_termination': ok_term,
        'total_energy_eh': float(m_e.group(1)) if m_e else None,
        'homo_lumo_gap_ev': float(m_g.group(1)) if m_g else None,
    }


def run_one(xyz: Path, outdir: Path, xtb_bin: str, timeout_s: int, gfn_level: int = 2):
    name = xyz.stem
    mdir = outdir / name
    mdir.mkdir(parents=True, exist_ok=True)

    inp = mdir / f"{name}.xyz"
    if not inp.exists():
        inp.write_text(Path(xyz).read_text(errors='ignore'))

    out = mdir / 'xtb.out'
    # NOTE:
    # - macOS xtb binary in this environment crashes in geometry optimization (--opt)
    #   due a Fortran runtime formatting error.
    # - Use single-point mode for robust pre-screening on existing 3D geometries.
    cmd = [xtb_bin, inp.name, '--gfn', str(gfn_level), '--sp']

    try:
        with out.open('w') as f:
            p = subprocess.run(cmd, cwd=mdir, stdout=f, stderr=subprocess.STDOUT, timeout=timeout_s)
        parsed = parse_xtb_out(out)
        ok = (p.returncode == 0) and parsed['normal_termination']
        detail = 'ok' if ok else f"returncode={p.returncode}"
        return ok, detail, parsed
    except subprocess.TimeoutExpired:
        return False, 'timeout', {'normal_termination': False, 'total_energy_eh': None, 'homo_lumo_gap_ev': None}
    except Exception as e:
        return False, f'error:{type(e).__name__}', {'normal_termination': False, 'total_energy_eh': None, 'homo_lumo_gap_ev': None}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--manifest', required=True, help='CSV with fields: idx,name,xyz_path')
    ap.add_argument('--batch-dir', required=True)
    ap.add_argument('--xtb-bin', default='/opt/homebrew/bin/xtb')
    ap.add_argument('--max-items', type=int, default=5000)
    ap.add_argument('--timeout', type=int, default=120)
    ap.add_argument('--gfn', type=int, default=2, choices=[0, 1, 2], help='xTB Hamiltonian level for --sp screening')
    args = ap.parse_args()

    manifest = Path(args.manifest)
    bdir = Path(args.batch_dir)
    outdir = bdir / 'xtb_results'
    outdir.mkdir(parents=True, exist_ok=True)

    state_file = bdir / 'xtb_state.json'
    result_csv = bdir / 'xtb_progress.csv'

    rows = list(csv.DictReader(manifest.open()))[:args.max_items]

    if state_file.exists():
        state = json.loads(state_file.read_text())
    else:
        state = {'next_index': 0, 'done': 0, 'ok': 0, 'fail': 0}

    start = state['next_index']
    print(f"start_index={start}, total={len(rows)}")

    write_header = not result_csv.exists()
    with result_csv.open('a', newline='') as rf:
        w = csv.writer(rf)
        if write_header:
            w.writerow([
                'idx', 'name', 'xyz_path', 'status', 'detail',
                'total_energy_eh', 'homo_lumo_gap_ev', 'normal_termination'
            ])

        for i in range(start, len(rows)):
            r = rows[i]
            xyz_path = Path(r['xyz_path'])
            ok, detail, parsed = run_one(xyz_path, outdir, args.xtb_bin, args.timeout, args.gfn)
            status = 'ok' if ok else 'fail'

            w.writerow([
                i + 1,
                r['name'],
                str(xyz_path),
                status,
                detail,
                parsed['total_energy_eh'],
                parsed['homo_lumo_gap_ev'],
                parsed['normal_termination'],
            ])
            rf.flush()

            state['done'] += 1
            if ok:
                state['ok'] += 1
            else:
                state['fail'] += 1
            state['next_index'] = i + 1
            state_file.write_text(json.dumps(state, indent=2))

            if (i + 1) % 50 == 0:
                print(f"progress {i+1}/{len(rows)} ok={state['ok']} fail={state['fail']}")

    print(f"finished done={state['done']} ok={state['ok']} fail={state['fail']}")


if __name__ == '__main__':
    main()
