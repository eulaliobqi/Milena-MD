#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Mapa de contatos receptor-ligante por resíduo, a partir da trajetória
GROMACS já processada (pós-equilíbrio).

O intervalo de resíduos do ligante é detectado da mesma forma que o módulo
ANALYSES (leitura direta da cadeia do complexo.pdb original) -- os
arquivos .gro/.tpr do GROMACS não preservam chain ID, só numeração de
resíduo, então a seleção usa resid, não chainID.

Uso:
  contact_map.py --complexo-pdb complexo.pdb --tpr md.tpr --xtc stable.xtc \
      --out-dir . [--ligand-chain B] [--cutoff 4.0]
"""
import argparse
import os
import sys

import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis import distances
from scipy import sparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def detect_ligand_range(complexo_pdb, ligand_chain='B'):
    first = last = None
    with open(complexo_pdb) as fh:
        for line in fh:
            if line.startswith(('ATOM', 'HETATM')) and line[21:22].strip() == ligand_chain:
                resnum = int(line[22:26])
                if first is None:
                    first = resnum
                last = resnum
    if first is None:
        raise SystemExit(f"ERRO: cadeia {ligand_chain} nao encontrada em {complexo_pdb}")
    return first, last


def residue_onehot(atomgroup):
    """Retorna (resids ordenados, matriz esparsa n_atoms x n_res one-hot)."""
    resids = sorted(set(int(r) for r in atomgroup.resids))
    idx = {r: i for i, r in enumerate(resids)}
    rows = np.arange(len(atomgroup))
    cols = np.array([idx[int(r)] for r in atomgroup.resids])
    data = np.ones(len(atomgroup), dtype=np.int32)
    onehot = sparse.csr_matrix((data, (rows, cols)), shape=(len(atomgroup), len(resids)))
    return resids, onehot


def resname_by_resid(atomgroup, resids):
    out = {}
    for r in resids:
        sel = atomgroup.select_atoms(f"resid {r}")
        out[r] = sel.resnames[0] if len(sel) else "?"
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--complexo-pdb", required=True)
    ap.add_argument("--tpr", required=True)
    ap.add_argument("--xtc", required=True)
    ap.add_argument("--ligand-chain", default="B")
    ap.add_argument("--cutoff", type=float, default=4.0,
                     help="Angstrom (default 4.0 = 0.4 nm, mesmo cutoff de ANALYSES)")
    ap.add_argument("--out-dir", default=".")
    args = ap.parse_args()

    lig_first, lig_last = detect_ligand_range(args.complexo_pdb, args.ligand_chain)
    print(f"[CONTACT_MAP] Ligante: residuos {lig_first}-{lig_last}", file=sys.stderr)

    u = mda.Universe(args.tpr, args.xtc)
    lig = u.select_atoms(f"resid {lig_first}-{lig_last} and not name H*")
    rec = u.select_atoms(f"protein and not (resid {lig_first}-{lig_last}) and not name H*")

    if len(lig) == 0 or len(rec) == 0:
        raise SystemExit(f"ERRO: selecao vazia (lig={len(lig)} atomos, rec={len(rec)} atomos)")

    rec_resids, rec_onehot = residue_onehot(rec)
    lig_resids, lig_onehot = residue_onehot(lig)
    rec_resname = resname_by_resid(rec, rec_resids)
    lig_resname = resname_by_resid(lig, lig_resids)

    n_res_rec, n_res_lig = len(rec_resids), len(lig_resids)
    contact_counts = np.zeros((n_res_rec, n_res_lig), dtype=np.int64)

    n_frames = len(u.trajectory)
    for ts in u.trajectory:
        d = distances.distance_array(rec.positions, lig.positions, box=ts.dimensions)
        atom_contact = sparse.csr_matrix(d < args.cutoff)
        # reduz contato átomo-átomo para "qualquer contato" por par de resíduos
        res_contact = (rec_onehot.T @ atom_contact @ lig_onehot) > 0
        contact_counts += res_contact.toarray()

    freq = contact_counts / max(n_frames, 1)

    os.makedirs(args.out_dir, exist_ok=True)

    # ── contact_map.csv ────────────────────────────────────────────────────
    rec_labels = [f"{rec_resname[r]}{r}" for r in rec_resids]
    lig_labels = [f"{lig_resname[r]}{r}" for r in lig_resids]
    csv_path = os.path.join(args.out_dir, "contact_map.csv")
    with open(csv_path, "w") as fh:
        fh.write("residue," + ",".join(lig_labels) + "\n")
        for i, rl in enumerate(rec_labels):
            fh.write(rl + "," + ",".join(f"{v:.4f}" for v in freq[i]) + "\n")
    print(f"[OK] {csv_path}", file=sys.stderr)

    # ── contact_map.png ────────────────────────────────────────────────────
    fig_h = max(6, n_res_rec * 0.12)
    fig_w = max(8, n_res_lig * 0.12)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    im = ax.imshow(freq, aspect='auto', cmap='viridis', vmin=0, vmax=1,
                    interpolation='nearest')
    ax.set_xlabel("Ligand residue")
    ax.set_ylabel("Receptor residue")
    ax.set_title("Receptor–Ligand Contact Frequency Map\n(fraction of frames, cutoff "
                  f"{args.cutoff/10:.1f} nm)")
    step_x = max(1, n_res_lig // 25)
    step_y = max(1, n_res_rec // 40)
    ax.set_xticks(range(0, n_res_lig, step_x))
    ax.set_xticklabels([lig_labels[i] for i in range(0, n_res_lig, step_x)],
                        rotation=90, fontsize=6)
    ax.set_yticks(range(0, n_res_rec, step_y))
    ax.set_yticklabels([rec_labels[i] for i in range(0, n_res_rec, step_y)], fontsize=6)
    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label("Contact frequency")
    plt.tight_layout()
    png_path = os.path.join(args.out_dir, "contact_map.png")
    plt.savefig(png_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"[OK] {png_path}", file=sys.stderr)

    # ── interface_residues.csv (para pharmacophore_profile.py) ────────────
    rec_max = freq.max(axis=1) if n_res_lig else np.zeros(n_res_rec)
    lig_max = freq.max(axis=0) if n_res_rec else np.zeros(n_res_lig)
    iface_path = os.path.join(args.out_dir, "interface_residues.csv")
    with open(iface_path, "w") as fh:
        fh.write("chain,resid,resname,max_contact_freq\n")
        for r, m in zip(rec_resids, rec_max):
            fh.write(f"receptor,{r},{rec_resname[r]},{m:.4f}\n")
        for r, m in zip(lig_resids, lig_max):
            fh.write(f"ligand,{r},{lig_resname[r]},{m:.4f}\n")
    print(f"[OK] {iface_path}", file=sys.stderr)
    print("Done.", file=sys.stderr)


if __name__ == "__main__":
    main()
