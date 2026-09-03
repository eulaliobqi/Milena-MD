#!/usr/bin/env python3
"""Lista residuos do receptor (cadeia A) com algum atomo pesado a <=4.5A de
qualquer atomo do peptideo (cadeia B), numa unica estrutura estatica
(pre-MD). Equivalente simplificado, para uma pose so, do que
bin/contact_map.py faz sobre trajetoria de MD (frequencia de contato).

Uso: interface_contacts.py <pdb> <out_csv> [--cutoff 4.5]
"""
import argparse
import csv
from pathlib import Path

from Bio.PDB import PDBParser, NeighborSearch


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("pdb")
    ap.add_argument("out_csv")
    ap.add_argument("--cutoff", type=float, default=4.5)
    args = ap.parse_args()

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("x", args.pdb)
    model = structure[0]

    chain_a_atoms = [a for a in model["A"].get_atoms() if a.element != "H"]
    chain_b_atoms = [a for a in model["B"].get_atoms() if a.element != "H"]

    ns = NeighborSearch(chain_a_atoms)
    contacted = {}
    for atom_b in chain_b_atoms:
        close = ns.search(atom_b.coord, args.cutoff)
        for atom_a in close:
            res = atom_a.get_parent()
            key = (res.id[1], res.resname)
            d = atom_a - atom_b
            if key not in contacted or d < contacted[key]:
                contacted[key] = d

    rows = sorted(({"resnum": k[0], "resname": k[1], "min_dist_A": round(v, 2)}
                    for k, v in contacted.items()), key=lambda r: r["resnum"])

    with open(args.out_csv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["resnum", "resname", "min_dist_A"])
        w.writeheader()
        w.writerows(rows)

    print(f"{len(rows)} residuos do receptor em contato (<= {args.cutoff} A) -> {args.out_csv}")
    print(", ".join(f"{r['resname']}{r['resnum']}" for r in rows))


if __name__ == "__main__":
    main()
