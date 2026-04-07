#!/usr/bin/env python3
"""
Build D/A topology-based candidate library with strict RDKit assembly rules.

Supported topologies:
- D-A
- D-A-D
- A-D-A
- D-pi-A
- D_n-A

This script is designed for the inquiry-driven workflow:
- topology selection is provided by initializer output
- initial sample count defaults to 10000
"""

from __future__ import annotations

import argparse
import csv
import itertools
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, List

from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")


@dataclass
class Frag:
    name: str
    smiles: str


def count_halogens(smiles: str) -> int:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return 0
    return sum(1 for a in mol.GetAtoms() if a.GetSymbol() in {"Cl", "Br", "I"})


def _copy_with_halogen_as_dummy(smiles: str, mapnums: List[int]) -> str:
    """Replace first N halogens with [*:mapnum] in order.

    If fragment already contains dummy mapped atoms, keep as-is when enough maps exist.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return ""

    existing = [a for a in mol.GetAtoms() if a.GetAtomicNum() == 0 and a.GetAtomMapNum() > 0]
    if len(existing) >= len(mapnums):
        return Chem.MolToSmiles(mol)

    rw = Chem.RWMol(mol)
    halogen_idx = [a.GetIdx() for a in rw.GetAtoms() if a.GetSymbol() in {"Cl", "Br", "I"}]
    if len(halogen_idx) < len(mapnums):
        return ""

    for idx, m in zip(halogen_idx[: len(mapnums)], mapnums):
        a = rw.GetAtomWithIdx(idx)
        a.SetAtomicNum(0)
        a.SetAtomMapNum(int(m))
        a.SetIsAromatic(False)

    out = rw.GetMol()
    Chem.SanitizeMol(out)
    return Chem.MolToSmiles(out)


def _molzip_join(parts: List[str]) -> str:
    combo = ".".join(parts)
    mol = Chem.MolFromSmiles(combo)
    if mol is None:
        return ""
    try:
        prod = Chem.molzip(mol)
    except Exception:
        return ""
    if prod is None:
        return ""
    smi = Chem.MolToSmiles(prod, canonical=True)
    if not smi or "." in smi:
        return ""
    return smi if Chem.MolFromSmiles(smi) is not None else ""


def assemble_d_a(d: Frag, a: Frag) -> str:
    d1 = _copy_with_halogen_as_dummy(d.smiles, [1])
    a1 = _copy_with_halogen_as_dummy(a.smiles, [1])
    if not d1 or not a1:
        return ""
    return _molzip_join([d1, a1])


def assemble_d_a_d(d1f: Frag, a: Frag, d2f: Frag) -> str:
    if count_halogens(a.smiles) < 2:
        return ""
    d1 = _copy_with_halogen_as_dummy(d1f.smiles, [1])
    d2 = _copy_with_halogen_as_dummy(d2f.smiles, [2])
    a12 = _copy_with_halogen_as_dummy(a.smiles, [1, 2])
    if not d1 or not d2 or not a12:
        return ""
    return _molzip_join([d1, a12, d2])


def assemble_a_d_a(a1f: Frag, d: Frag, a2f: Frag) -> str:
    if count_halogens(d.smiles) < 2:
        return ""
    a1 = _copy_with_halogen_as_dummy(a1f.smiles, [1])
    a2 = _copy_with_halogen_as_dummy(a2f.smiles, [2])
    d12 = _copy_with_halogen_as_dummy(d.smiles, [1, 2])
    if not a1 or not a2 or not d12:
        return ""
    return _molzip_join([a1, d12, a2])


def assemble_d_pi_a(d: Frag, a: Frag, pi_smiles: str) -> str:
    # pi bridge should provide two mapped positions
    d1 = _copy_with_halogen_as_dummy(d.smiles, [1])
    a2 = _copy_with_halogen_as_dummy(a.smiles, [2])
    p12 = _copy_with_halogen_as_dummy(pi_smiles, [1, 2])
    if not d1 or not a2 or not p12:
        return ""
    return _molzip_join([d1, p12, a2])


def assemble_d_n_a(donors: List[Frag], a: Frag) -> str:
    n = len(donors)
    if n < 3:
        return ""
    if count_halogens(a.smiles) < n:
        return ""
    parts = []
    for i, d in enumerate(donors, 1):
        di = _copy_with_halogen_as_dummy(d.smiles, [i])
        if not di:
            return ""
        parts.append(di)
    ai = _copy_with_halogen_as_dummy(a.smiles, list(range(1, n + 1)))
    if not ai:
        return ""
    parts.append(ai)
    return _molzip_join(parts)


def read_frags(path: Path) -> List[Frag]:
    rows = []
    if path.suffix.lower() == ".smi":
        for i, line in enumerate(path.read_text().splitlines(), 1):
            s = line.strip()
            if s:
                rows.append(Frag(name=f"frag_{i}", smiles=s.split()[0]))
        return rows

    with path.open() as f:
        r = csv.DictReader(f)
        for i, row in enumerate(r, 1):
            smi = (row.get("smiles") or "").strip()
            if not smi:
                continue
            nm = (row.get("name") or f"frag_{i}").strip()
            rows.append(Frag(name=nm, smiles=smi))
    return rows


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--donors", required=True, help="Path to donor .csv/.smi")
    ap.add_argument("--acceptors", required=True, help="Path to acceptor .csv/.smi")
    ap.add_argument("--topology", required=True, choices=["D-A", "D-A-D", "A-D-A", "D-pi-A", "D_n-A"])
    ap.add_argument("--initial-sample-count", type=int, default=10000)
    ap.add_argument("--pi-bridge-smiles", default="Brc1ccccc1Br")
    ap.add_argument("--dn", type=int, default=3, help="n for D_n-A topology")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    donors = read_frags(Path(args.donors))
    acceptors = read_frags(Path(args.acceptors))

    out_rows = []

    if args.topology == "D-A":
        for d, a in itertools.product(donors, acceptors):
            smi = assemble_d_a(d, a)
            if smi:
                out_rows.append([d.name, a.name, args.topology, smi])
            if len(out_rows) >= args.initial_sample_count:
                break

    elif args.topology == "D-A-D":
        for d1, a, d2 in itertools.product(donors, acceptors, donors):
            smi = assemble_d_a_d(d1, a, d2)
            if smi:
                out_rows.append([f"{d1.name}|{d2.name}", a.name, args.topology, smi])
            if len(out_rows) >= args.initial_sample_count:
                break

    elif args.topology == "A-D-A":
        for a1, d, a2 in itertools.product(acceptors, donors, acceptors):
            smi = assemble_a_d_a(a1, d, a2)
            if smi:
                out_rows.append([d.name, f"{a1.name}|{a2.name}", args.topology, smi])
            if len(out_rows) >= args.initial_sample_count:
                break

    elif args.topology == "D-pi-A":
        for d, a in itertools.product(donors, acceptors):
            smi = assemble_d_pi_a(d, a, args.pi_bridge_smiles)
            if smi:
                out_rows.append([d.name, a.name, args.topology, smi])
            if len(out_rows) >= args.initial_sample_count:
                break

    else:  # D_n-A
        n = max(3, int(args.dn))
        for a in acceptors:
            for donor_tuple in itertools.combinations_with_replacement(donors, n):
                smi = assemble_d_n_a(list(donor_tuple), a)
                if smi:
                    dnames = "|".join(d.name for d in donor_tuple)
                    out_rows.append([dnames, a.name, f"D_{n}-A", smi])
                if len(out_rows) >= args.initial_sample_count:
                    break
            if len(out_rows) >= args.initial_sample_count:
                break

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["donor", "acceptor", "topology", "product_smiles"])
        w.writerows(out_rows)

    print(f"generated={len(out_rows)} topology={args.topology} out={out}")


if __name__ == "__main__":
    main()
