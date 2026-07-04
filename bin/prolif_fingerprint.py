#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Fingerprint de interação (ProLIF) receptor-ligante ao longo da trajetória.

Módulo de maior risco de dependência do conjunto (prolif+rdkit+MDAnalysis é
a combinação menos testada neste projeto) -- por isso escreve um log de
diagnóstico ANTES de tentar rodar o ProLIF e cai num fallback com CSV vazio
em caso de erro, em vez de estourar traceback cru (mesmo padrão defensivo
já usado em MMGBSA_ROBUST).

Uso:
  prolif_fingerprint.py --complexo-pdb complexo.pdb --tpr md.tpr \
      --xtc stable.xtc --out-dir . [--ligand-chain B] [--stride 1]
"""
import argparse
import os
import sys
import traceback


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


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--complexo-pdb", required=True)
    ap.add_argument("--tpr", required=True)
    ap.add_argument("--xtc", required=True)
    ap.add_argument("--ligand-chain", default="B")
    ap.add_argument("--stride", type=int, default=1,
                     help="Usa 1 a cada N frames (default 1 = todos)")
    ap.add_argument("--out-dir", default=".")
    args = ap.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    log_path = os.path.join(args.out_dir, "prolif_log.txt")
    csv_path = os.path.join(args.out_dir, "prolif_fingerprint.csv")
    png_path = os.path.join(args.out_dir, "prolif_heatmap.png")

    def log(msg):
        print(msg, file=sys.stderr)
        with open(log_path, "a") as fh:
            fh.write(msg + "\n")

    log("=== PROLIF_FINGERPRINT ===")

    try:
        import MDAnalysis as mda
        import prolif as plf
        import pandas as pd
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        lig_first, lig_last = detect_ligand_range(args.complexo_pdb, args.ligand_chain)
        log(f"[OK] Ligante: residuos {lig_first}-{lig_last}")

        u = mda.Universe(args.tpr, args.xtc)
        lig_sel = u.select_atoms(f"resid {lig_first}-{lig_last}")
        rec_sel = u.select_atoms(f"protein and not (resid {lig_first}-{lig_last})")
        log(f"[OK] Selecoes: ligante={len(lig_sel)} atomos, receptor={len(rec_sel)} atomos")

        fp = plf.Fingerprint([
            "HBAcceptor", "HBDonor", "Hydrophobic", "PiStacking",
            "Anionic", "Cationic", "VdWContact",
        ])
        traj = u.trajectory[::args.stride]
        log(f"[INFO] Rodando ProLIF sobre {len(traj)} frames (stride={args.stride})...")
        fp.run(traj, lig_sel, rec_sel)

        df = fp.to_dataframe()
        df.to_csv(csv_path)
        log(f"[OK] {csv_path} ({df.shape[0]} frames x {df.shape[1]} interacoes)")

        # Frequência por (residuo, tipo_interacao) ao longo da trajetória
        freq = df.mean(axis=0)
        freq_df = freq.reset_index()
        freq_df.columns = ["ligand_residue", "protein_residue", "interaction", "frequency"]
        pivot = freq_df.pivot_table(index="protein_residue", columns="interaction",
                                     values="frequency", aggfunc="sum", fill_value=0.0)

        if pivot.empty:
            raise ValueError("Nenhuma interacao detectada pelo ProLIF (pivot vazio)")

        fig, ax = plt.subplots(figsize=(8, max(4, len(pivot) * 0.25)))
        im = ax.imshow(pivot.values, aspect='auto', cmap='magma', vmin=0, vmax=1)
        ax.set_xticks(range(len(pivot.columns)))
        ax.set_xticklabels(pivot.columns, rotation=45, ha='right', fontsize=8)
        ax.set_yticks(range(len(pivot.index)))
        ax.set_yticklabels(pivot.index, fontsize=7)
        ax.set_title("ProLIF: Frequência de Interação por Resíduo do Receptor")
        fig.colorbar(im, ax=ax, label="Frequência")
        plt.tight_layout()
        plt.savefig(png_path, dpi=150, bbox_inches='tight')
        plt.close()
        log(f"[OK] {png_path}")
        log("[MMGBSA-STYLE OK] ProLIF concluido com sucesso")

    except Exception:
        log("ERRO: ProLIF falhou -- ver traceback abaixo")
        log(traceback.format_exc())
        # Fallback -- nao bloqueia downstream (mesmo padrao de MMGBSA_ROBUST)
        if not os.path.exists(csv_path):
            with open(csv_path, "w") as fh:
                fh.write("ligand_residue,protein_residue,interaction,frequency\n")
        log("[FALLBACK] CSV vazio gravado -- ProLIF nao disponivel/falhou nesta rodada")


if __name__ == "__main__":
    main()
