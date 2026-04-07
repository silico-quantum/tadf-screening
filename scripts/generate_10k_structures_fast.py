#!/usr/bin/env python3
from __future__ import annotations

import csv
import json
import shutil
from pathlib import Path
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem

RDLogger.DisableLog('rdApp.*')

ROOT = Path('/Users/molbot/.openclaw/workspace/tadf-screening')
VALID_PRODUCTS = ROOT / 'datasets' / 'blue_library_valid_products.csv'
PLAN_CSV = ROOT / 'datasets' / 'blue_screening_plan_10k.csv'
POOL_DIR = ROOT / 'structures' / 'blue_pool_xyz'
PLAN_DIR = ROOT / 'structures' / 'blue_plan_10k_xyz'
MANIFEST = ROOT / 'datasets' / 'blue_plan_10k_manifest.csv'
SUMMARY = ROOT / 'datasets' / 'blue_plan_10k_summary.json'
TARGET = 10000


def safe_name(sample_id: str) -> str:
    return sample_id.replace('::', '__').replace('/', '_')


def load_smiles_map() -> dict[str, str]:
    rows = list(csv.DictReader(VALID_PRODUCTS.open()))
    return {r['name']: r['product_smiles'] for r in rows if r.get('status') == 'ok' and r.get('product_smiles')}


def smiles_to_xyz_quick(smiles: str) -> str:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return ''
    mol = Chem.AddHs(mol)

    params = AllChem.ETKDGv3()
    params.randomSeed = 0xF00D
    params.maxIterations = 500

    rc = AllChem.EmbedMolecule(mol, params)
    if rc != 0:
        rc = AllChem.EmbedMolecule(mol, useRandomCoords=True, randomSeed=0xF00D)
        if rc != 0:
            return ''

    try:
        AllChem.UFFOptimizeMolecule(mol, maxIters=200)
    except Exception:
        pass

    try:
        return Chem.MolToXYZBlock(mol)
    except Exception:
        return ''


def count_existing() -> int:
    return sum(1 for _ in PLAN_DIR.glob('*.xyz'))


def main():
    POOL_DIR.mkdir(parents=True, exist_ok=True)
    PLAN_DIR.mkdir(parents=True, exist_ok=True)

    smiles_map = load_smiles_map()
    plans = list(csv.DictReader(PLAN_CSV.open()))

    done = count_existing()
    missing_base = 0
    xyz_failed = 0
    reused_base = 0
    built_base = 0

    with MANIFEST.open('w', newline='') as mf:
        w = csv.writer(mf)
        w.writerow(['idx', 'sample_id', 'name', 'status', 'detail', 'xyz_path'])

        for idx, row in enumerate(plans, start=1):
            if done >= TARGET:
                break

            sample_id = row['sample_id']
            name = row['name']
            out_xyz = PLAN_DIR / f"{safe_name(sample_id)}.xyz"

            if out_xyz.exists() and out_xyz.stat().st_size > 0:
                done += 1
                w.writerow([idx, sample_id, name, 'ok', 'exists', str(out_xyz)])
                continue

            smi = smiles_map.get(name)
            if not smi:
                missing_base += 1
                w.writerow([idx, sample_id, name, 'skip', 'no_valid_product_smiles', ''])
                continue

            base_xyz = POOL_DIR / f"{name}.xyz"
            if not (base_xyz.exists() and base_xyz.stat().st_size > 0):
                xyz = smiles_to_xyz_quick(smi)
                if not xyz:
                    xyz_failed += 1
                    w.writerow([idx, sample_id, name, 'skip', 'smiles_to_xyz_failed', ''])
                    continue
                base_xyz.write_text(xyz)
                built_base += 1
            else:
                reused_base += 1

            shutil.copy2(base_xyz, out_xyz)
            done += 1
            w.writerow([idx, sample_id, name, 'ok', 'copied_from_base', str(out_xyz)])

            if done % 250 == 0:
                print(f'progress {done}/{TARGET}', flush=True)

    summary = {
        'target': TARGET,
        'done': done,
        'missing_base': missing_base,
        'xyz_failed': xyz_failed,
        'built_base': built_base,
        'reused_base': reused_base,
    }
    SUMMARY.write_text(json.dumps(summary, ensure_ascii=False, indent=2))
    print(json.dumps(summary, ensure_ascii=False), flush=True)


if __name__ == '__main__':
    main()
