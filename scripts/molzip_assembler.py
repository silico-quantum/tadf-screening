#!/usr/bin/env python3
"""RDKit molzip-based fragment assembler.

Inspired by the R-group/mapped-fragment workflow:
- input fragments carry atom-map labels (e.g., [*:1], [*:2])
- fragments are joined via Chem.molzip after dot-joining

Examples:
  python tools/molzip_assembler.py one \
    --frags "[*:1]c1ccccc1" "[*:1]N1c2ccccc2Oc2ccccc21"

  python tools/molzip_assembler.py build \
    --core "[*:1]c1nc([*:2])nc([*:3])n1" \
    --r1 data/r1.txt --r2 data/r2.txt --r3 data/r3.txt \
    --out datasets/molzip_products.csv
"""

from __future__ import annotations

import argparse
import csv
import itertools
from pathlib import Path
from typing import Iterable

from rdkit import Chem, RDLogger

RDLogger.DisableLog('rdApp.*')


def sanitize_smiles(smi: str) -> str:
    return smi.strip()


def read_smiles_file(path: Path) -> list[str]:
    rows = []
    for line in path.read_text().splitlines():
        s = sanitize_smiles(line)
        if s:
            rows.append(s)
    return rows


def try_molzip(frag_smiles: list[str]) -> tuple[bool, str]:
    """Return (ok, product_or_reason)."""
    joined = '.'.join(frag_smiles)
    mol = Chem.MolFromSmiles(joined)
    if mol is None:
        return False, 'invalid_fragment_combo'

    try:
        prod = Chem.molzip(mol)
    except Exception as e:
        return False, f'molzip_failed:{type(e).__name__}'

    if prod is None:
        return False, 'molzip_none'

    try:
        smi = Chem.MolToSmiles(prod, canonical=True)
    except Exception as e:
        return False, f'to_smiles_failed:{type(e).__name__}'

    if not smi or '.' in smi:
        return False, 'disconnected_or_empty'

    if Chem.MolFromSmiles(smi) is None:
        return False, 'invalid_product_smiles'

    return True, smi


def cmd_one(args):
    ok, out = try_molzip(args.frags)
    print('ok' if ok else 'fail', out)


def cmd_build(args):
    core = args.core
    r1 = read_smiles_file(Path(args.r1))
    r2 = read_smiles_file(Path(args.r2))
    r3 = read_smiles_file(Path(args.r3)) if args.r3 else ['']

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    rows = []
    total = 0
    ok_n = 0

    for a, b, c in itertools.product(r1, r2, r3):
        total += 1
        frags = [core, a, b] + ([c] if c else [])
        frags = [x for x in frags if x]
        ok, prod = try_molzip(frags)
        row = {
            'idx': total,
            'status': 'ok' if ok else 'fail',
            'core': core,
            'r1': a,
            'r2': b,
            'r3': c,
            'product_smiles': prod if ok else '',
            'reason': '' if ok else prod,
        }
        if ok:
            ok_n += 1
        rows.append(row)

    with out_path.open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    print(f'total={total} ok={ok_n} fail={total-ok_n} out={out_path}')


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest='cmd', required=True)

    p1 = sub.add_parser('one')
    p1.add_argument('--frags', nargs='+', required=True)

    p2 = sub.add_parser('build')
    p2.add_argument('--core', required=True)
    p2.add_argument('--r1', required=True)
    p2.add_argument('--r2', required=True)
    p2.add_argument('--r3')
    p2.add_argument('--out', required=True)

    args = ap.parse_args()
    if args.cmd == 'one':
        cmd_one(args)
    else:
        cmd_build(args)


if __name__ == '__main__':
    main()
