#!/usr/bin/env python3
"""Reproduz os numeros da secao 7 de comparativo_vs_DN2954-GORE12T.md
(divergencia de conformacao HADDOCK3 vs Boltz-2 para os peptideos novos).

Tres checagens, sobre arquivos ja commitados em
novas-analises/{DN2954-GORE1-2T,DN2954-GORE1-2T-GGS3}/haddock/ e
novas-analises/diagnostico_divergencia_pose/:

1. RMSD (Kabsch, backbone N/CA/C/O da cadeia B) da pose HADDOCK original
   ("best_pose_haddock.pdb", run1-guided_orig_frozen-peptide) e da pose
   HADDOCK reprocessada com flexref cobrindo o peptideo inteiro
   ("*_flexfull_best.pdb", run1-guided) contra os 3 conformeros de
   entrada do ensemble (helice/fita/PPII) -- mostra que os dois runs
   HADDOCK mantem o backbone essencialmente identico ao conformero
   helice de entrada.
2. Distancia ponta-a-ponta (Ca-Ca) e raio de giro de cada pose, incluindo
   as 3 amostras de difusao do Boltz-2 -- mostra convergencia (GGS3) ou
   divergencia (GORE1-2T) entre amostras independentes do mesmo modelo.
3. pLDDT medio por amostra Boltz-2 (B-factor da cadeia B no PDB).

Uso: python scripts/diagnose_pose_divergence.py
"""
from pathlib import Path

import numpy as np
from Bio.PDB import PDBParser

BASE = Path(__file__).resolve().parents[1]
DIAG = BASE / "diagnostico_divergencia_pose"
P = PDBParser(QUIET=True)


def bb_coords(path, chain="B"):
    s = P.get_structure("x", path)
    ch = s[0][chain]
    coords = []
    for res in ch:
        if res.id[0] != " ":
            continue
        for name in ("N", "CA", "C", "O"):
            if name in res:
                coords.append(res[name].coord)
    return np.array(coords)


def ca_stats(path, chain="B"):
    s = P.get_structure("x", path)
    ch = s[0][chain]
    cas, bfacs = [], []
    for res in ch:
        if res.id[0] == " " and "CA" in res:
            cas.append(res["CA"].coord)
            bfacs.append(res["CA"].get_bfactor())
    coords = np.array(cas)
    end_to_end = np.linalg.norm(coords[0] - coords[-1])
    rg = np.sqrt(np.mean(np.sum((coords - coords.mean(axis=0)) ** 2, axis=1)))
    return len(cas), end_to_end, rg, np.mean(bfacs)


def kabsch_rmsd(P_, Q_):
    Pc, Qc = P_ - P_.mean(axis=0), Q_ - Q_.mean(axis=0)
    H = Pc.T @ Qc
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    D = np.diag([1, 1, d])
    R = Vt.T @ D @ U.T
    return float(np.sqrt(np.mean(np.sum(((R @ Pc.T).T - Qc) ** 2, axis=1))))


SYSTEMS = {
    "gore12t": {
        "haddock_orig": BASE / "DN2954-GORE1-2T/haddock/best_pose_haddock.pdb",
        "haddock_flexfull": DIAG / "haddock_gore12t_flexfull_best.pdb",
        "helix": DIAG / "gore1-2t_helix.pdb",
        "strand": DIAG / "gore1-2t_strand.pdb",
        "ppii": DIAG / "gore1-2t_ppii.pdb",
        "boltz_models": [DIAG / f"dn2954_gore1-2t_model_{i}.pdb" for i in range(3)],
    },
    "ggs3": {
        "haddock_orig": BASE / "DN2954-GORE1-2T-GGS3/haddock/best_pose_haddock.pdb",
        "haddock_flexfull": DIAG / "haddock_ggs3_flexfull_best.pdb",
        "helix": DIAG / "gore1-2t-ggs3_helix.pdb",
        "strand": DIAG / "gore1-2t-ggs3_strand.pdb",
        "ppii": DIAG / "gore1-2t-ggs3_ppii.pdb",
        "boltz_models": [DIAG / f"dn2954_gore1-2t-ggs3_model_{i}.pdb" for i in range(3)],
    },
}

for tag, paths in SYSTEMS.items():
    print(f"\n=== {tag} ===")
    bb_orig = bb_coords(paths["haddock_orig"])
    bb_full = bb_coords(paths["haddock_flexfull"])
    for conf_name in ("helix", "strand", "ppii"):
        bb_ref = bb_coords(paths[conf_name])
        n = min(len(bb_orig), len(bb_ref))
        r_orig = kabsch_rmsd(bb_orig[:n], bb_ref[:n])
        r_full = kabsch_rmsd(bb_full[:n], bb_ref[:n])
        print(f"  RMSD vs input {conf_name:8s}: HADDOCK original={r_orig:5.2f} A | "
              f"HADDOCK flexref-full={r_full:5.2f} A")

    for label, path in (("HADDOCK original", paths["haddock_orig"]),
                         ("HADDOCK flexref-full", paths["haddock_flexfull"])):
        n, e2e, rg, _ = ca_stats(path)
        print(f"  {label:22s}: n={n} end-to-end={e2e:6.1f} A  Rg={rg:5.1f} A")

    for i, mpath in enumerate(paths["boltz_models"]):
        n, e2e, rg, plddt = ca_stats(mpath)
        print(f"  Boltz-2 model_{i}       : n={n} end-to-end={e2e:6.1f} A  "
              f"Rg={rg:5.1f} A  pLDDT_medio={plddt:5.1f}")
