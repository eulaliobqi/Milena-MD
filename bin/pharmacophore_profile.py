#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Perfil farmacofórico dos resíduos de interface (persistentes ao longo da
trajetória), a partir de interface_residues.csv gerado por contact_map.py.

Não usa RDKit: é uma interface peptídeo-peptídeo (proteína-proteína), não
molécula pequena -- classificação por tipo de resíduo (tabela fixa) é
suficiente e muito mais robusta que qualquer perception química genérica.

Uso:
  pharmacophore_profile.py --interface-csv interface_residues.csv \
      --out-dir . [--threshold 0.2]
"""
import argparse
import os
import sys
from collections import Counter

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Normaliza variantes de protonação/nomenclatura (AMBER/CHARMM) para o
# código de 3 letras canônico.
_RESNAME_ALIAS = {
    "HSD": "HIS", "HSE": "HIS", "HSP": "HIS", "HID": "HIS", "HIE": "HIS", "HIP": "HIS",
    "ASH": "ASP", "GLH": "GLU", "LYN": "LYS", "CYX": "CYS", "CYM": "CYS",
}

# Features farmacofóricas por resíduo -- tabela determinística.
_FEATURES = {
    "ALA": {"hydrophobic"},
    "ARG": {"positive_ionizable", "hbond_donor"},
    "ASN": {"hbond_donor", "hbond_acceptor"},
    "ASP": {"negative_ionizable", "hbond_acceptor"},
    "CYS": {"hydrophobic"},
    "GLN": {"hbond_donor", "hbond_acceptor"},
    "GLU": {"negative_ionizable", "hbond_acceptor"},
    "GLY": set(),  # sem cadeia lateral
    "HIS": {"aromatic", "hbond_donor", "hbond_acceptor", "positive_ionizable"},
    "ILE": {"hydrophobic"},
    "LEU": {"hydrophobic"},
    "LYS": {"positive_ionizable", "hbond_donor"},
    "MET": {"hydrophobic"},
    "PHE": {"aromatic", "hydrophobic"},
    "PRO": {"hydrophobic"},
    "SER": {"hbond_donor", "hbond_acceptor"},
    "THR": {"hbond_donor", "hbond_acceptor"},
    "TRP": {"aromatic", "hydrophobic", "hbond_donor"},
    "TYR": {"aromatic", "hydrophobic", "hbond_donor", "hbond_acceptor"},
    "VAL": {"hydrophobic"},
}

_FEATURE_ORDER = ["hydrophobic", "aromatic", "hbond_donor", "hbond_acceptor",
                   "positive_ionizable", "negative_ionizable"]


def canonical(resname):
    resname = resname.upper()
    return _RESNAME_ALIAS.get(resname, resname)


def classify(resname):
    return _FEATURES.get(canonical(resname), set())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--interface-csv", required=True)
    ap.add_argument("--out-dir", default=".")
    ap.add_argument("--threshold", type=float, default=0.2,
                     help="Frequência mínima de contato para considerar 'persistente' (default 0.2)")
    args = ap.parse_args()

    rows = []
    with open(args.interface_csv) as fh:
        header = fh.readline().strip().split(",")
        for ln in fh:
            parts = ln.strip().split(",")
            if len(parts) != len(header):
                continue
            rec = dict(zip(header, parts))
            rec["max_contact_freq"] = float(rec["max_contact_freq"])
            rows.append(rec)

    persistent = [r for r in rows if r["max_contact_freq"] >= args.threshold]
    print(f"[PHARMACOPHORE] {len(persistent)}/{len(rows)} residuos acima do "
          f"limiar de frequencia {args.threshold}", file=sys.stderr)

    os.makedirs(args.out_dir, exist_ok=True)

    # ── pharmacophore_profile.csv ──────────────────────────────────────────
    csv_path = os.path.join(args.out_dir, "pharmacophore_profile.csv")
    with open(csv_path, "w") as fh:
        fh.write("chain,resid,resname,max_contact_freq,features\n")
        for r in sorted(persistent, key=lambda x: -x["max_contact_freq"]):
            feats = classify(r["resname"])
            fh.write(f"{r['chain']},{r['resid']},{r['resname']},"
                     f"{r['max_contact_freq']:.4f},\"{';'.join(sorted(feats)) or 'none'}\"\n")
    print(f"[OK] {csv_path}", file=sys.stderr)

    # ── pharmacophore_summary.png ──────────────────────────────────────────
    counts_by_chain = {"receptor": Counter(), "ligand": Counter()}
    for r in persistent:
        for feat in classify(r["resname"]):
            counts_by_chain[r["chain"]][feat] += 1

    x = range(len(_FEATURE_ORDER))
    w = 0.35
    fig, ax = plt.subplots(figsize=(9, 5))
    rec_vals = [counts_by_chain["receptor"].get(f, 0) for f in _FEATURE_ORDER]
    lig_vals = [counts_by_chain["ligand"].get(f, 0) for f in _FEATURE_ORDER]
    ax.bar([i - w / 2 for i in x], rec_vals, width=w, label="Receptor",
           color="steelblue", edgecolor="black", linewidth=0.5)
    ax.bar([i + w / 2 for i in x], lig_vals, width=w, label="Ligand",
           color="tomato", edgecolor="black", linewidth=0.5)
    ax.set_xticks(list(x))
    ax.set_xticklabels([f.replace("_", " ") for f in _FEATURE_ORDER], rotation=20, ha='right')
    ax.set_ylabel("N. residuos de interface")
    ax.set_title(f"Perfil Farmacofórico da Interface\n(resíduos com contato ≥ {args.threshold:.0%} dos frames)")
    ax.legend()
    ax.grid(alpha=0.25, axis='y')
    plt.tight_layout()
    png_path = os.path.join(args.out_dir, "pharmacophore_summary.png")
    plt.savefig(png_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"[OK] {png_path}", file=sys.stderr)
    print("Done.", file=sys.stderr)


if __name__ == "__main__":
    main()
