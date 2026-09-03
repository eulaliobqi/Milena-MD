#!/usr/bin/env python3
"""Protona um complexo (2 cadeias, sem HETATM nao-padrao) em pH 8,2 via
pdb2pqr30 --titration-state-method propka, e extrai do PROPKA especificamente
HIS/ASP/GLU/CYS com pKa previsto, sinalizando anomalias em relacao ao pH 8,2.

Reaproveita a receita e o bug-fix ja documentados em
eulalio-pos-doc/codigo/fase9_blocoF6b_protonation/{protonate_receptors,protonate_complexes}.py:
a tabela de pKa do PROPKA sai em STDERR (nao STDOUT) e so' deve ser lida do
combined stdout+stderr.

Uso: protonate_ph82.py <input.pdb> <out_dir> <label> [--pH 8.2]
Saidas em <out_dir>/: <label>_ph8.2.pdb, <label>_propka.log, <label>_pka_resumo.csv
"""
import argparse
import csv
import re
import subprocess
import sys
from pathlib import Path

TITULAVEIS = {"HIS", "ASP", "GLU", "CYS", "HID", "HIE", "HIP", "ASH", "GLH", "CYX"}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("input_pdb")
    ap.add_argument("out_dir")
    ap.add_argument("label")
    ap.add_argument("--pH", type=float, default=8.2)
    args = ap.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    out_pdb = out_dir / f"{args.label}_ph{args.pH}.pdb"
    out_pqr = out_dir / f"{args.label}_ph{args.pH}.pqr"
    log_path = out_dir / f"{args.label}_propka.log"
    csv_path = out_dir / f"{args.label}_pka_resumo.csv"

    cmd = [
        "pdb2pqr30", "--ff", "AMBER", "--titration-state-method", "propka",
        "--with-ph", str(args.pH), "--keep-chain", "--pdb-output", str(out_pdb),
        args.input_pdb, str(out_pqr),
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    combined = result.stdout + "\n---STDERR---\n" + result.stderr
    log_path.write_text(combined)

    ok = result.returncode == 0 and out_pdb.exists()
    if not ok:
        sys.exit(f"ERRO: pdb2pqr30 falhou para {args.input_pdb} (rc={result.returncode}). Ver {log_path}")

    # extrai tabela de pKa (formato "SUMMARY OF THIS PREDICTION")
    rows = []
    if "SUMMARY OF THIS PREDICTION" in combined:
        block = combined.split("SUMMARY OF THIS PREDICTION")[1]
        for m in re.finditer(r"^\s*([A-Z]{2,4})\s+(\d+)\s+([A-Za-z])\s+([\d.]+)\s+([\d.]+)",
                              block, re.MULTILINE):
            resname, resnum, chain, pka, model_pka = m.groups()
            if resname not in ("HIS", "ASP", "GLU", "CYS"):
                continue
            pka = float(pka)
            anomalia = ""
            if resname == "HIS":
                estado = "protonada/carregada (HIP) - INCOMUM em pH 8,2" if pka > args.pH else "neutra (HID/HIE)"
                if pka > args.pH:
                    anomalia = "pKa acima do pH - confirmar tautomero/protonacao"
            elif resname == "CYS":
                estado = "possivel tiolato (desprotonada)" if pka < args.pH else "SH neutro"
                if abs(pka - args.pH) < 0.3:
                    anomalia = "pKa limitrofe (~pH 8,2) - estado incerto, revisar manualmente"
            else:  # ASP/GLU
                estado = "carregada (padrao)" if pka < args.pH else "protonada (ASH/GLH) - INCOMUM em pH 8,2"
                if pka > args.pH:
                    anomalia = "pKa acima do pH - efeito local anomalo, revisar"
            rows.append({
                "resname": resname, "resnum": resnum, "chain": chain,
                "pKa_previsto": pka, "pH_referencia": args.pH,
                "estado_atribuido": estado, "anomalia": anomalia,
            })

    with open(csv_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["resname", "resnum", "chain", "pKa_previsto",
                                           "pH_referencia", "estado_atribuido", "anomalia"])
        w.writeheader()
        w.writerows(rows)

    n_anom = sum(1 for r in rows if r["anomalia"])
    print(f"OK: {args.label} protonado em pH {args.pH} -> {out_pdb}")
    print(f"    {len(rows)} residuos titulaveis (HIS/ASP/GLU/CYS) reportados, {n_anom} com anomalia.")
    print(f"    Resumo: {csv_path}")


if __name__ == "__main__":
    main()
