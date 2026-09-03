#!/usr/bin/env python3
"""Roda PRODIGY-prot (Xue et al. 2016) sobre uma lista de complexos
protonados (cadeia A=receptor, B=peptideo), extraindo dG e Kd previstos.

Ressalva conhecida, ja documentada no grupo (ver
eulalio-pos-doc/codigo/fase9_blocoF6c_docking_energies/run_prodigy.py):
PRODIGY foi calibrado em interfaces proteina-proteina "normais" (dominios
estruturados); validade estatistica para peptideo curto/nao-globular nao
esta comprovada - reportado do mesmo jeito, mas marcado na coluna `caveat`.

Uso: run_prodigy_batch.py <out_csv> <label1>=<pdb1> [<label2>=<pdb2> ...]
"""
import csv
import re
import subprocess
import sys

CAVEAT = "PRODIGY calibrado em proteina-proteina estruturada; peptideo curto/nao-globular fora do regime validado"


def main():
    if len(sys.argv) < 3:
        sys.exit("Uso: run_prodigy_batch.py <out_csv> label=pdb [label=pdb ...]")
    out_csv = sys.argv[1]
    items = [a.split("=", 1) for a in sys.argv[2:]]

    rows = []
    for label, pdb_path in items:
        try:
            result = subprocess.run(["prodigy", pdb_path, "--selection", "A", "B"],
                                     capture_output=True, text=True, timeout=180)
            m_dg = re.search(r"Predicted binding affinity \(kcal\.mol-1\):\s*(-?[\d.]+)", result.stdout)
            m_kd = re.search(r"Predicted dissociation constant \(M\).*:\s*([\d.eE+-]+)", result.stdout)
            dg = float(m_dg.group(1)) if m_dg else None
            kd = float(m_kd.group(1)) if m_kd else None
            status = "OK" if dg is not None else f"SEM_MATCH_rc={result.returncode}"
        except Exception as exc:
            dg, kd, status = None, None, f"EXCECAO: {exc}"
        rows.append({"label": label, "pdb": pdb_path, "delta_g_kcal_mol": dg,
                     "kd_M_25C": kd, "status": status, "caveat": CAVEAT})
        print(f"{label}: dG={dg} kcal/mol, Kd={kd} M - {status}")

    with open(out_csv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["label", "pdb", "delta_g_kcal_mol", "kd_M_25C", "status", "caveat"])
        w.writeheader()
        w.writerows(rows)

    n_ok = sum(1 for r in rows if r["status"] == "OK")
    print(f"\n{n_ok}/{len(rows)} energias PRODIGY calculadas com sucesso -> {out_csv}")


if __name__ == "__main__":
    main()
