#!/usr/bin/env python3
"""Robust SMILES assembly utilities for TADF screening.

Goals:
1) Validate donor/acceptor fragment SMILES before assembly
2) Assemble connected products (prefer D-A, fallback D-A-D)
3) Reject malformed/disconnected products
4) Provide CLI for audit + build steps used by screening workflow

Usage examples:
  python tools/smiles_assembler.py audit \
      --donors data/donors_expanded.json --acceptors data/acceptors_expanded.json \
      --out datasets/fragment_audit.csv

  python tools/smiles_assembler.py build \
      --candidates datasets/blue_library_candidates_full.csv \
      --donors data/donors_expanded.json --acceptors data/acceptors_expanded.json \
      --out datasets/blue_library_valid_products.csv
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Dict, Tuple, Optional

from rdkit import Chem, RDLogger

RDLogger.DisableLog('rdApp.*')


def is_valid_smiles(smiles: str) -> bool:
    if not smiles or not isinstance(smiles, str):
        return False
    try:
        return Chem.MolFromSmiles(smiles) is not None
    except Exception:
        return False


def canonical_smiles(smiles: str) -> Optional[str]:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return Chem.MolToSmiles(mol, canonical=True)


def has_donor_attach_site(smiles: str) -> bool:
    """Heuristic: donor should have N (preferably ring N with attachable environment)."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 7:  # N
            return True
    return False


def has_acceptor_attach_site(smiles: str) -> bool:
    """Heuristic: acceptor should have aromatic C attach point."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6 and atom.GetIsAromatic():
            return True
    return False


def load_fragments(path: Path, expected_type: str) -> Dict[str, dict]:
    data = json.loads(path.read_text())
    out: Dict[str, dict] = {}
    for f in data.get('fragments', []):
        if f.get('type') != expected_type:
            continue
        abbr = f.get('abbr')
        if not abbr:
            continue
        out[abbr] = f
    return out


def _load_connect_fragments(project_root: Path):
    import importlib.util

    mb_path = project_root / 'tools' / 'molecule_builder.py'
    spec = importlib.util.spec_from_file_location('molecule_builder', mb_path)
    mod = importlib.util.module_from_spec(spec)
    assert spec and spec.loader
    spec.loader.exec_module(mod)
    return mod.connect_fragments


def assemble_connected_smiles(
    donor_smiles: str,
    acceptor_smiles: str,
    connect_fragments,
) -> Tuple[Optional[str], str]:
    """Try to assemble connected product.

    Strategy:
    - Try D-A first (closer to current candidate definition)
    - Fallback D-A-D if D-A fails
    - Reject disconnected product (contains '.')
    """
    for pattern in ('D-A', 'D-A-D'):
        try:
            mol = connect_fragments(donor_smiles, acceptor_smiles, pattern=pattern)
        except Exception:
            mol = None

        if mol is None:
            continue

        try:
            smi = Chem.MolToSmiles(mol, canonical=True)
        except Exception:
            continue

        if not smi or '.' in smi:
            continue

        if Chem.MolFromSmiles(smi) is None:
            continue

        return smi, f'ok:{pattern}'

    return None, 'assembly_failed'


def cmd_audit(args):
    donors = load_fragments(Path(args.donors), 'donor')
    acceptors = load_fragments(Path(args.acceptors), 'acceptor')

    rows = []

    for abbr, d in donors.items():
        smi = d.get('smiles', '')
        rows.append({
            'kind': 'donor',
            'abbr': abbr,
            'name': d.get('name', ''),
            'smiles': smi,
            'valid_smiles': is_valid_smiles(smi),
            'attach_ok': has_donor_attach_site(smi) if is_valid_smiles(smi) else False,
            'canonical_smiles': canonical_smiles(smi) if is_valid_smiles(smi) else '',
        })

    for abbr, a in acceptors.items():
        smi = a.get('smiles', '')
        rows.append({
            'kind': 'acceptor',
            'abbr': abbr,
            'name': a.get('name', ''),
            'smiles': smi,
            'valid_smiles': is_valid_smiles(smi),
            'attach_ok': has_acceptor_attach_site(smi) if is_valid_smiles(smi) else False,
            'canonical_smiles': canonical_smiles(smi) if is_valid_smiles(smi) else '',
        })

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    total = len(rows)
    valid = sum(1 for r in rows if r['valid_smiles'])
    attach_ok = sum(1 for r in rows if r['attach_ok'])
    print(f'audit_total={total} valid={valid} attach_ok={attach_ok} out={out}')


def cmd_build(args):
    project_root = Path(args.project_root)
    connect_fragments = _load_connect_fragments(project_root)

    donors = load_fragments(Path(args.donors), 'donor')
    acceptors = load_fragments(Path(args.acceptors), 'acceptor')

    cands = list(csv.DictReader(Path(args.candidates).open()))
    rows = []

    for r in cands:
        name = r['name']
        d_abbr = r['donor']
        a_abbr = r['acceptor']
        d = donors.get(d_abbr)
        a = acceptors.get(a_abbr)

        if d is None or a is None:
            rows.append({
                'name': name,
                'donor': d_abbr,
                'acceptor': a_abbr,
                'status': 'missing_fragment',
                'product_smiles': '',
                'reason': 'donor_or_acceptor_not_found',
            })
            continue

        d_smi = d.get('smiles', '')
        a_smi = a.get('smiles', '')

        if not is_valid_smiles(d_smi):
            rows.append({
                'name': name,
                'donor': d_abbr,
                'acceptor': a_abbr,
                'status': 'invalid',
                'product_smiles': '',
                'reason': 'invalid_donor_smiles',
            })
            continue

        if not is_valid_smiles(a_smi):
            rows.append({
                'name': name,
                'donor': d_abbr,
                'acceptor': a_abbr,
                'status': 'invalid',
                'product_smiles': '',
                'reason': 'invalid_acceptor_smiles',
            })
            continue

        product, msg = assemble_connected_smiles(d_smi, a_smi, connect_fragments)
        if product is None:
            rows.append({
                'name': name,
                'donor': d_abbr,
                'acceptor': a_abbr,
                'status': 'invalid',
                'product_smiles': '',
                'reason': msg,
            })
            continue

        rows.append({
            'name': name,
            'donor': d_abbr,
            'acceptor': a_abbr,
            'status': 'ok',
            'product_smiles': product,
            'reason': msg,
        })

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    total = len(rows)
    ok = sum(1 for r in rows if r['status'] == 'ok')
    fail = total - ok
    print(f'build_total={total} ok={ok} fail={fail} out={out}')


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest='cmd', required=True)

    ap_audit = sub.add_parser('audit')
    ap_audit.add_argument('--donors', required=True)
    ap_audit.add_argument('--acceptors', required=True)
    ap_audit.add_argument('--out', required=True)

    ap_build = sub.add_parser('build')
    ap_build.add_argument('--project-root', default='/Users/molbot/.openclaw/workspace/tadf-screening')
    ap_build.add_argument('--candidates', required=True)
    ap_build.add_argument('--donors', required=True)
    ap_build.add_argument('--acceptors', required=True)
    ap_build.add_argument('--out', required=True)

    args = ap.parse_args()
    if args.cmd == 'audit':
        cmd_audit(args)
    else:
        cmd_build(args)


if __name__ == '__main__':
    main()
