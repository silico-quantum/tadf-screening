#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import math
import re
import shutil
import subprocess
from pathlib import Path

EV_RE = re.compile(r"([-+]?\d*\.\d+|\d+)\s*eV", re.IGNORECASE)


def first_excitation_ev_from_text(text: str):
    # Heuristic parser: first positive value followed by 'eV'
    vals = []
    for m in EV_RE.finditer(text):
        try:
            v = float(m.group(1))
            if v > 0:
                vals.append(v)
        except Exception:
            pass
    return vals[0] if vals else None


def run_shell(cmd_tmpl: str, cwd: Path, xyz_name: str, timeout: int):
    cmd = cmd_tmpl.format(xyz=xyz_name)
    p = subprocess.run(
        cmd,
        cwd=str(cwd),
        shell=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        timeout=timeout,
    )
    return p.returncode, p.stdout


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--xtb-progress", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--target-nm", type=float, default=470.0)
    ap.add_argument("--window-nm", type=float, default=30.0)
    ap.add_argument("--max-items", type=int, default=500)
    ap.add_argument("--timeout", type=int, default=180)
    ap.add_argument("--prep-cmd", default="xtb4stda {xyz}", help="prep command template")
    ap.add_argument("--run-cmd", default="stda -xtb -e 10", help="TDDFT-xTB command template")
    args = ap.parse_args()

    xtb_progress = Path(args.xtb_progress)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    csv_out = outdir / "tddft_xtb_results.csv"
    pass_out = outdir / "tddft_xtb_blue_window.csv"

    rows = list(csv.DictReader(xtb_progress.open()))
    ok_rows = [r for r in rows if r.get("status") == "ok"][: args.max_items]

    # Tool availability check
    prep_tool = args.prep_cmd.split()[0]
    run_tool = args.run_cmd.split()[0]
    prep_ok = shutil.which(prep_tool) is not None
    run_ok = shutil.which(run_tool) is not None

    results = []
    for i, r in enumerate(ok_rows, 1):
        name = r.get("name", f"mol_{i}")
        xyz_path = Path(r["xyz_path"])
        wdir = outdir / "work" / name
        wdir.mkdir(parents=True, exist_ok=True)

        local_xyz = wdir / f"{name}.xyz"
        if not local_xyz.exists():
            local_xyz.write_text(xyz_path.read_text(errors="ignore"))

        if not (prep_ok and run_ok):
            results.append(
                {
                    "idx": i,
                    "name": name,
                    "xyz_path": str(xyz_path),
                    "status": "fail",
                    "detail": f"tool_missing: prep={prep_tool}:{prep_ok}, run={run_tool}:{run_ok}",
                    "s1_ev": "",
                    "lambda_nm": "",
                    "target_nm": args.target_nm,
                    "delta_nm": "",
                    "pass_blue_window": "",
                }
            )
            continue

        try:
            rc1, out1 = run_shell(args.prep_cmd, wdir, local_xyz.name, args.timeout)
            (wdir / "prep.log").write_text(out1)
            if rc1 != 0:
                results.append(
                    {
                        "idx": i,
                        "name": name,
                        "xyz_path": str(xyz_path),
                        "status": "fail",
                        "detail": f"prep_failed:{rc1}",
                        "s1_ev": "",
                        "lambda_nm": "",
                        "target_nm": args.target_nm,
                        "delta_nm": "",
                        "pass_blue_window": "",
                    }
                )
                continue

            rc2, out2 = run_shell(args.run_cmd, wdir, local_xyz.name, args.timeout)
            (wdir / "stda.log").write_text(out2)
            if rc2 != 0:
                results.append(
                    {
                        "idx": i,
                        "name": name,
                        "xyz_path": str(xyz_path),
                        "status": "fail",
                        "detail": f"tddft_failed:{rc2}",
                        "s1_ev": "",
                        "lambda_nm": "",
                        "target_nm": args.target_nm,
                        "delta_nm": "",
                        "pass_blue_window": "",
                    }
                )
                continue

            s1_ev = first_excitation_ev_from_text(out2)
            if s1_ev is None or s1_ev <= 0:
                results.append(
                    {
                        "idx": i,
                        "name": name,
                        "xyz_path": str(xyz_path),
                        "status": "fail",
                        "detail": "no_excitation_parsed",
                        "s1_ev": "",
                        "lambda_nm": "",
                        "target_nm": args.target_nm,
                        "delta_nm": "",
                        "pass_blue_window": "",
                    }
                )
                continue

            lam = 1239.84 / s1_ev
            delta = abs(lam - args.target_nm)
            keep = delta <= args.window_nm

            results.append(
                {
                    "idx": i,
                    "name": name,
                    "xyz_path": str(xyz_path),
                    "status": "ok",
                    "detail": "ok",
                    "s1_ev": round(s1_ev, 6),
                    "lambda_nm": round(lam, 2),
                    "target_nm": args.target_nm,
                    "delta_nm": round(delta, 2),
                    "pass_blue_window": keep,
                }
            )

        except subprocess.TimeoutExpired:
            results.append(
                {
                    "idx": i,
                    "name": name,
                    "xyz_path": str(xyz_path),
                    "status": "fail",
                    "detail": "timeout",
                    "s1_ev": "",
                    "lambda_nm": "",
                    "target_nm": args.target_nm,
                    "delta_nm": "",
                    "pass_blue_window": "",
                }
            )

    fields = [
        "idx",
        "name",
        "xyz_path",
        "status",
        "detail",
        "s1_ev",
        "lambda_nm",
        "target_nm",
        "delta_nm",
        "pass_blue_window",
    ]

    with csv_out.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(results)

    passed = [r for r in results if str(r.get("pass_blue_window")).lower() == "true"]
    with pass_out.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(passed)

    n_ok = sum(1 for r in results if r["status"] == "ok")
    n_fail = len(results) - n_ok
    print(
        f"total={len(results)} ok={n_ok} fail={n_fail} blue_pass={len(passed)} out={csv_out}"
    )


if __name__ == "__main__":
    main()
