#!/usr/bin/env python3
"""Calcula, estaticamente (sem MD), as distancias-chave da triade catalitica
His79/Asp126/Ser222 + Asp216 (S1, numeracao ja usada em
assets/samplesheet_dn2954_gore12t.csv) numa pose de complexo (receptor cadeia
A), para triagem de qualidade de pose antes de prosseguir para
protonacao/interacao/afinidade (ANALYSES_TRIAD do Milena-MD faz o mesmo
sobre trajetoria de MD; aqui e' a mesma logica de distancia aplicada a uma
unica estrutura estatica).

Uso: triad_distances.py <pdb> <label> [--chain A] [--his 79] [--asp 126] [--ser 222] [--asp1 216]
"""
import argparse
import csv
import sys
from pathlib import Path

try:
    from Bio.PDB import PDBParser
except ImportError:
    sys.exit("ERRO: biopython nao disponivel neste ambiente (precisa de Bio.PDB).")


def get_atom(structure, chain_id, resnum, atom_name):
    for model in structure:
        for chain in model:
            if chain.id != chain_id:
                continue
            for res in chain:
                if res.id[1] == resnum:
                    if atom_name in res:
                        return res[atom_name]
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("pdb")
    ap.add_argument("label")
    ap.add_argument("--chain", default="A")
    ap.add_argument("--his", type=int, default=79)
    ap.add_argument("--asp", type=int, default=126)
    ap.add_argument("--ser", type=int, default=222)
    ap.add_argument("--asp1", type=int, default=216)
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("x", args.pdb)

    pairs = [
        ("His-NE2...Ser-OG", args.his, "NE2", args.ser, "OG"),
        ("His-ND1...Asp-OD1", args.his, "ND1", args.asp, "OD1"),
        ("His-ND1...Asp-OD2", args.his, "ND1", args.asp, "OD2"),
        ("Ser-OG...Asp216-OD1", args.ser, "OG", args.asp1, "OD1"),
        ("Ser-OG...Asp216-OD2", args.ser, "OG", args.asp1, "OD2"),
    ]

    rows = []
    for name, r1, a1, r2, a2 in pairs:
        atom1 = get_atom(structure, args.chain, r1, a1)
        atom2 = get_atom(structure, args.chain, r2, a2)
        if atom1 is None or atom2 is None:
            rows.append({"par": name, "distancia_A": None, "status": "ATOMO_AUSENTE"})
            continue
        d = atom1 - atom2
        status = "OK (<=4.5A)" if d <= 4.5 else "DISTANTE (>4.5A, revisar pose)"
        rows.append({"par": name, "distancia_A": round(float(d), 2), "status": status})

    out_path = Path(args.out) if args.out else Path(args.pdb).with_suffix("").with_name(f"{args.label}_triad_distances.csv")
    with open(out_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["par", "distancia_A", "status"])
        w.writeheader()
        w.writerows(rows)

    print(f"{args.label}:")
    for r in rows:
        print(f"  {r['par']}: {r['distancia_A']} A - {r['status']}")
    print(f"  -> {out_path}")


if __name__ == "__main__":
    main()
