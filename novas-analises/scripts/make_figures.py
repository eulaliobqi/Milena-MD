#!/usr/bin/env python3
"""Gera as figuras do Bloco 1 (DN2954 x GORE1-2T / GORE1-2T-GGS3) a partir
dos CSVs/JSONs ja produzidos pelas etapas anteriores. Salva em novas-analises/figuras/.
"""
import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = Path("/home/eulalio/gromacs/Milena-MD/novas-analises")
OUT = BASE / "figuras"
OUT.mkdir(exist_ok=True)

SYS_LABEL = {"DN2954-GORE1-2T": "GORE1-2T (21 aa)", "DN2954-GORE1-2T-GGS3": "GORE1-2T(GGS)3 (75 aa)"}
COLORS = {"DN2954-GORE1-2T": "#1f77b4", "DN2954-GORE1-2T-GGS3": "#d62728"}

# ---------- Fig 1: HADDOCK cluster scores ----------
fig, axes = plt.subplots(1, 2, figsize=(10, 4.2), sharey=False)
for ax, sysname in zip(axes, SYS_LABEL):
    tsv = BASE / sysname / "haddock" / "run1-guided" / "10_caprieval" / "capri_clt.tsv"
    rows = []
    with open(tsv) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.split("\t")
            if parts[0] == "cluster_rank":
                header = parts
                continue
            rows.append(dict(zip(header, parts)))
    ranks = [int(r["cluster_rank"]) for r in rows]
    scores = [float(r["score"]) for r in rows]
    stds = [float(r["score_std"]) for r in rows]
    ax.bar(ranks, scores, yerr=stds, color=COLORS[sysname], alpha=0.85, capsize=3)
    ax.set_title(SYS_LABEL[sysname], fontsize=10)
    ax.set_xlabel("cluster (rank)")
    ax.set_ylabel("HADDOCK score")
    ax.axhline(0, color="grey", lw=0.6)
fig.suptitle("Fig. 1 — HADDOCK3 (padrão-ouro): score por cluster, 09_seletopclusts", fontsize=11)
fig.tight_layout()
fig.savefig(OUT / "fig1_haddock_scores.png", dpi=160)
plt.close(fig)

# ---------- Fig 2: Boltz-2 confidence per sample ----------
conf_paths = {
    "DN2954-GORE1-2T": BASE / "DN2954-GORE1-2T/boltz2/boltz_out/boltz_results_dn2954_gore1-2t/predictions/dn2954_gore1-2t",
    "DN2954-GORE1-2T-GGS3": BASE / "DN2954-GORE1-2T-GGS3/boltz2/boltz_out/boltz_results_dn2954_gore1-2t-ggs3/predictions/dn2954_gore1-2t-ggs3",
}
metrics = ["confidence_score", "iptm", "ptm", "complex_plddt"]
fig, ax = plt.subplots(figsize=(9, 4.5))
width = 0.1
x0 = 0
xticks, xlabels = [], []
for sysname, d in conf_paths.items():
    for m_idx in range(3):
        f = d / f"confidence_dn2954_{'gore1-2t' if sysname=='DN2954-GORE1-2T' else 'gore1-2t-ggs3'}_model_{m_idx}.json"
        data = json.load(open(f))
        for j, metric in enumerate(metrics):
            ax.bar(x0 + j * width, data[metric], width=width,
                   color=plt.cm.tab10(j), label=metric if x0 == 0 else None)
        xticks.append(x0 + width * 1.5)
        xlabels.append(f"{SYS_LABEL[sysname].split(' ')[0]}\nmodel_{m_idx}")
        x0 += width * 4 + 0.15
ax.set_xticks(xticks)
ax.set_xticklabels(xlabels, fontsize=8)
ax.set_ylabel("valor (0-1)")
ax.set_ylim(0, 1.05)
ax.legend(fontsize=8, ncol=4, loc="lower center", bbox_to_anchor=(0.5, -0.28))
fig.suptitle("Fig. 2 — Boltz-2 (fronteira): métricas de confiança por amostra de difusão", fontsize=11)
fig.tight_layout()
fig.savefig(OUT / "fig2_boltz_confidence.png", dpi=160, bbox_inches="tight")
plt.close(fig)

# ---------- Fig 3: Triad distances consistency ----------
fig, ax = plt.subplots(figsize=(9, 4.5))
pose_labels, ne2_ser, nd1_asp = [], [], []
for sysname in SYS_LABEL:
    short = "gore1-2t" if sysname == "DN2954-GORE1-2T" else "gore1-2t-ggs3"
    haddock_csv = BASE / sysname / "qualidade" / "haddock_triad.csv"
    rows = list(csv.DictReader(open(haddock_csv)))
    d = {r["par"]: float(r["distancia_A"]) for r in rows}
    pose_labels.append(f"{short}\nHADDOCK")
    ne2_ser.append(d["His-NE2...Ser-OG"])
    nd1_asp.append(min(d["His-ND1...Asp-OD1"], d["His-ND1...Asp-OD2"]))
    for m in range(3):
        boltz_csv = BASE / sysname / "qualidade" / f"boltz_triad_model{m}.csv"
        rows = list(csv.DictReader(open(boltz_csv)))
        d = {r["par"]: float(r["distancia_A"]) for r in rows}
        pose_labels.append(f"{short}\nBoltz-2 m{m}")
        ne2_ser.append(d["His-NE2...Ser-OG"])
        nd1_asp.append(min(d["His-ND1...Asp-OD1"], d["His-ND1...Asp-OD2"]))

x = range(len(pose_labels))
ax.plot(x, ne2_ser, "o-", label="His-NE2···Ser-OG", color="#2ca02c")
ax.plot(x, nd1_asp, "s-", label="His-ND1···Asp-OD (menor)", color="#9467bd")
ax.axhline(4.5, color="grey", ls="--", lw=1, label="limiar de aceitação (4,5 Å)")
ax.set_xticks(list(x))
ax.set_xticklabels(pose_labels, rotation=45, ha="right", fontsize=8)
ax.set_ylabel("distância (Å)")
ax.legend(fontsize=8)
fig.suptitle("Fig. 3 — Geometria da díade catalítica His79/Ser222/Asp126, todas as poses", fontsize=11)
fig.tight_layout()
fig.savefig(OUT / "fig3_triad_distances.png", dpi=160)
plt.close(fig)

# ---------- Fig 4: pKa summary (HIS/ASP/GLU), highlight anomaly ----------
fig, axes = plt.subplots(1, 4, figsize=(13, 4.2), sharey=True)
pose_files = [
    ("GORE1-2T\nHADDOCK", BASE / "DN2954-GORE1-2T/protonacao_ph8.2/haddock_gore1-2t_pka_resumo.csv"),
    ("GORE1-2T\nBoltz-2 m0", BASE / "DN2954-GORE1-2T/protonacao_ph8.2/boltz_gore1-2t_model0_pka_resumo.csv"),
    ("GGS3\nHADDOCK", BASE / "DN2954-GORE1-2T-GGS3/protonacao_ph8.2/haddock_gore1-2t-ggs3_pka_resumo.csv"),
    ("GGS3\nBoltz-2 m0", BASE / "DN2954-GORE1-2T-GGS3/protonacao_ph8.2/boltz_gore1-2t-ggs3_model0_pka_resumo.csv"),
]
color_by_res = {"HIS": "#1f77b4", "ASP": "#d62728", "GLU": "#ff7f0e"}
for ax, (label, path) in zip(axes, pose_files):
    rows = [r for r in csv.DictReader(open(path)) if r["resname"] != "CYS"]
    for r in rows:
        pka = float(r["pKa_previsto"])
        c = color_by_res[r["resname"]]
        marker = "*" if r["anomalia"] else "o"
        size = 140 if r["anomalia"] else 40
        ax.scatter(r["resname"], pka, color=c, marker=marker, s=size,
                   edgecolor="black" if r["anomalia"] else "none", zorder=3)
    ax.axhline(8.2, color="grey", ls="--", lw=1)
    ax.set_title(label, fontsize=9)
axes[0].set_ylabel("pKa previsto (PROPKA)")
fig.suptitle("Fig. 4 — pKa de HIS/ASP/GLU nas 4 poses protonadas (linha tracejada = pH 8,2; ★ = anomalia)", fontsize=10.5)
fig.tight_layout()
fig.savefig(OUT / "fig4_pka_summary.png", dpi=160)
plt.close(fig)

# ---------- Fig 5: PRODIGY affinity ----------
fig, ax = plt.subplots(figsize=(7, 4.5))
bars_x, bars_y, bars_c = [], [], []
labels = []
i = 0
for sysname in SYS_LABEL:
    csv_path = BASE / sysname / "afinidade" / "prodigy_results.csv"
    rows = list(csv.DictReader(open(csv_path)))
    for r in rows:
        bars_x.append(i)
        bars_y.append(float(r["delta_g_kcal_mol"]))
        bars_c.append(COLORS[sysname])
        labels.append(f"{SYS_LABEL[sysname].split(' ')[0]}\n{r['label']}")
        i += 1
ax.bar(bars_x, bars_y, color=bars_c)
ax.set_xticks(bars_x)
ax.set_xticklabels(labels, fontsize=8)
ax.set_ylabel("ΔG previsto (kcal/mol)")
ax.axhline(0, color="black", lw=0.8)
fig.suptitle("Fig. 5 — Afinidade prevista (PRODIGY-prot), poses protonadas pH 8,2", fontsize=11)
fig.tight_layout()
fig.savefig(OUT / "fig5_prodigy_affinity.png", dpi=160)
plt.close(fig)

# ---------- Fig 6: Interface overlap vs baseline MD ----------
baseline_top = {55,56,57,58,63,99,120,121,122,123,166,167,170,175,197,218,219,220,240,242}
fig, ax = plt.subplots(figsize=(7, 4.5))
bars_x, bars_y, bars_c, labels = [], [], [], []
i = 0
for sysname in SYS_LABEL:
    for method in ["haddock", "boltz"]:
        csv_path = BASE / sysname / "interacao" / f"{method}_interface.csv"
        rows = list(csv.DictReader(open(csv_path)))
        resnums = {int(r["resnum"]) for r in rows}
        overlap = len(resnums & baseline_top)
        bars_x.append(i)
        bars_y.append(overlap)
        bars_c.append(COLORS[sysname])
        labels.append(f"{SYS_LABEL[sysname].split(' ')[0]}\n{method}")
        i += 1
ax.bar(bars_x, bars_y, color=bars_c)
ax.axhline(len(baseline_top), color="grey", ls="--", lw=1, label=f"total no baseline MD (n={len(baseline_top)})")
ax.set_xticks(bars_x)
ax.set_xticklabels(labels, fontsize=8)
ax.set_ylabel("nº resíduos em comum com a MD de 100 ns (DN2954-GORE12T)")
ax.legend(fontsize=8)
fig.suptitle("Fig. 6 — Recuperação da interface já validada por MD, por método", fontsize=10.5)
fig.tight_layout()
fig.savefig(OUT / "fig6_interface_overlap.png", dpi=160)
plt.close(fig)

print("Figuras geradas em", OUT)
for p in sorted(OUT.glob("*.png")):
    print(" -", p.name)
