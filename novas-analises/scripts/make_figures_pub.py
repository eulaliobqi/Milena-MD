#!/usr/bin/env python3
"""Regenera as figuras do Bloco 1 em resolucao de publicacao (PDF vetorial +
PNG >=600dpi, estilo Nature de partida), usando os helpers do skill
scientific-visualization (style_presets.py, figure_export.py) em vez de
parametros de figura escolhidos a mao.

Escopo corrigido em 2026-09-03: mostra so o que segue para a proxima etapa
(pose Boltz-2, decisao registrada em artigo.md secao 2.2/2.7 e
comparativo_vs_DN2954-GORE12T.md secao 7). O HADDOCK3 continua descrito
nos Metodos como validacao cruzada do sitio/interface, mas seus valores
de trade/protonacao/afinidade/interface (poses de peptideo nao-refinadas,
ver secao 7 do comparativo) nao entram mais nas figuras principais -- os
numeros completos com HADDOCK ficam so no comparativo, nao aqui.

Roda localmente (Windows), lendo os dados ja sincronizados de volta do
servidor em novas-analises/. Nao afirma conformidade de submissao --
apenas aplica o snapshot de estilo/DPI do skill.
"""
import csv
import json
import sys
from pathlib import Path

SKILL_SCRIPTS = Path(r"C:\Users\eulal\.claude\skills\scientific-visualization\scripts")
sys.path.insert(0, str(SKILL_SCRIPTS))

from style_presets import style_context, figure_size_for_profile  # noqa: E402
from figure_export import export_figure  # noqa: E402

import matplotlib.pyplot as plt  # noqa: E402

BASE = Path(__file__).resolve().parents[1]
OUT = BASE / "figuras"
OUT.mkdir(exist_ok=True)

SYS_LABEL = {"DN2954-GORE1-2T": "GORE1-2T (21 aa)", "DN2954-GORE1-2T-GGS3": "GORE1-2T(GGS)3 (75 aa)"}
COLOR_KEY = {"DN2954-GORE1-2T": 0, "DN2954-GORE1-2T-GGS3": 1}

PROVENANCE_COMMON = {
    "raw_data": "novas-analises/{DN2954-GORE1-2T,DN2954-GORE1-2T-GGS3}/{qualidade,protonacao_ph8.2,interacao,afinidade,boltz2}/",
    "transformations": ["nenhuma alem da leitura direta dos CSV/TSV/JSON ja produzidos pelas etapas do Bloco 1",
                         "poses HADDOCK excluidas destas figuras por decisao de escopo (nao seguem para a MD; "
                         "ver comparativo_vs_DN2954-GORE12T.md secao 7 para os valores completos com HADDOCK)"],
}


def mm_figsize(width_mm, height_mm):
    return (width_mm / 25.4, height_mm / 25.4)


def export(fig, name, provenance_extra=None):
    prov = dict(PROVENANCE_COMMON)
    if provenance_extra:
        prov.update(provenance_extra)
    report = export_figure(
        fig, OUT / name, formats=["pdf", "png"], dpi=600,
        bbox_inches=None, font_mode="truetype", overwrite=True,
        write_manifest=True, provenance=prov,
    )
    plt.close(fig)
    for o in report["outputs"]:
        print(" ->", o["path"], o["size_bytes"], "bytes")


# ---------- Fig 1: Boltz-2 confidence per sample ----------
conf_paths = {
    "DN2954-GORE1-2T": (BASE / "DN2954-GORE1-2T/boltz2/boltz_out/boltz_results_dn2954_gore1-2t/predictions/dn2954_gore1-2t", "gore1-2t"),
    "DN2954-GORE1-2T-GGS3": (BASE / "DN2954-GORE1-2T-GGS3/boltz2/boltz_out/boltz_results_dn2954_gore1-2t-ggs3/predictions/dn2954_gore1-2t-ggs3", "gore1-2t-ggs3"),
}
metrics = ["confidence_score", "iptm", "ptm", "complex_plddt"]
metric_labels = {"confidence_score": "confidence", "iptm": "ipTM", "ptm": "pTM", "complex_plddt": "pLDDT"}
with style_context("nature", palette_name="okabe_ito_on_white") as info:
    print("style:", info)
    w, h = mm_figsize(136, 80)
    fig, ax = plt.subplots(figsize=(w, h), layout="constrained")
    width = 0.19
    x0 = 0
    xticks, xlabels = [], []
    palette = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    for sysname, (d, tag) in conf_paths.items():
        for m_idx in range(3):
            f = d / f"confidence_dn2954_{tag}_model_{m_idx}.json"
            data = json.load(open(f))
            for j, metric in enumerate(metrics):
                ax.bar(x0 + j * width, data[metric], width=width, color=palette[j % len(palette)],
                       edgecolor="black", linewidth=0.3,
                       label=metric_labels[metric] if x0 == 0 else None)
            xticks.append(x0 + width * 1.5)
            xlabels.append(f"{SYS_LABEL[sysname].split(' ')[0]}\nmodel_{m_idx}")
            x0 += width * 4 + 0.18
    ax.set_xticks(xticks)
    ax.set_xticklabels(xlabels)
    ax.set_ylabel("valor (0-1)")
    ax.set_ylim(0, 1.05)
    ax.legend(ncol=4, loc="upper center", bbox_to_anchor=(0.5, -0.22), frameon=False)
    export(fig, "fig1_boltz_confidence",
           {"raw_data_detail": "confidence_*.json de cada amostra de difusao (Boltz-2 --diffusion_samples 3)"})

# ---------- Fig 2: Triad distances consistency (Boltz-2 only, 3 samples x 2 systems) ----------
with style_context("nature", palette_name="okabe_ito_on_white") as info:
    w, h = mm_figsize(120, 78)
    fig, ax = plt.subplots(figsize=(w, h), layout="constrained")
    pose_labels, ne2_ser, nd1_asp = [], [], []
    for sysname in SYS_LABEL:
        short = "gore1-2t" if sysname == "DN2954-GORE1-2T" else "gore1-2t-ggs3"
        for m in range(3):
            boltz_csv = BASE / sysname / "qualidade" / f"boltz_triad_model{m}.csv"
            rows = list(csv.DictReader(open(boltz_csv)))
            d = {r["par"]: float(r["distancia_A"]) for r in rows}
            pose_labels.append(f"{short}\nBoltz-2 m{m}")
            ne2_ser.append(d["His-NE2...Ser-OG"])
            nd1_asp.append(min(d["His-ND1...Asp-OD1"], d["His-ND1...Asp-OD2"]))
    x = list(range(len(pose_labels)))
    palette = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    ax.plot(x, ne2_ser, marker="o", label="His-NE2···Ser-OG", color=palette[2])
    ax.plot(x, nd1_asp, marker="s", label="His-ND1···Asp-OD (menor)", color=palette[3])
    ax.axhline(4.5, color="0.4", ls="--", lw=0.7, label="limiar de aceitacao (4,5 A)")
    ax.set_xticks(x)
    ax.set_xticklabels(pose_labels, rotation=45, ha="right")
    ax.set_ylabel("distancia (A)")
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.42), ncol=3, frameon=False)
    export(fig, "fig2_triad_distances",
           {"uncertainty": "medida unica por pose (estrutura estatica, sem replicata de MD nesta etapa)"})

# ---------- Fig 3: pKa summary (HIS/ASP/GLU), Boltz-2 model_0 only, 2 panels ----------
with style_context("nature", palette_name="okabe_ito_on_white") as info:
    w, h = mm_figsize(100, 62)
    fig, axes = plt.subplots(1, 2, figsize=(w, h), sharey=True, layout="constrained")
    pose_files = [
        ("GORE1-2T\nBoltz-2 m0", BASE / "DN2954-GORE1-2T/protonacao_ph8.2/boltz_gore1-2t_model0_pka_resumo.csv"),
        ("GGS3\nBoltz-2 m0", BASE / "DN2954-GORE1-2T-GGS3/protonacao_ph8.2/boltz_gore1-2t-ggs3_model0_pka_resumo.csv"),
    ]
    palette = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    color_by_res = {"HIS": palette[0], "ASP": palette[1], "GLU": palette[4 % len(palette)]}
    for ax, (label, path), panel in zip(axes, pose_files, "ab"):
        rows = [r for r in csv.DictReader(open(path)) if r["resname"] != "CYS"]
        for r in rows:
            pka = float(r["pKa_previsto"])
            c = color_by_res[r["resname"]]
            anomaly = bool(r["anomalia"])
            ax.scatter(r["resname"], pka, color=c, marker="*" if anomaly else "o",
                       s=90 if anomaly else 22, edgecolor="black" if anomaly else "none",
                       linewidth=0.6, zorder=3)
        ax.axhline(8.2, color="0.4", ls="--", lw=0.7)
        ax.set_title(label)
        ax.text(-0.18, 1.06, panel, transform=ax.transAxes, fontsize=8,
                fontweight="bold", va="bottom", ha="right")
    axes[0].set_ylabel("pKa previsto (PROPKA)")
    export(fig, "fig3_pka_summary",
           {"missing_data": "CYS omitida do painel (pKa=99,99 sentinela de nao-titulavel/dissulfeto em todas as poses, ver texto)"})

# ---------- Fig 4: PRODIGY affinity (Boltz-2 only) ----------
with style_context("nature", palette_name="okabe_ito_on_white") as info:
    w, h = mm_figsize(75, 65)
    fig, ax = plt.subplots(figsize=(w, h), layout="constrained")
    palette = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    bars_x, bars_y, bars_c, labels = [], [], [], []
    i = 0
    for sysname in SYS_LABEL:
        rows = list(csv.DictReader(open(BASE / sysname / "afinidade" / "prodigy_results.csv")))
        for r in rows:
            if r["label"] == "haddock":
                continue
            bars_x.append(i)
            bars_y.append(float(r["delta_g_kcal_mol"]))
            bars_c.append(palette[COLOR_KEY[sysname]])
            labels.append(SYS_LABEL[sysname].split(" ")[0])
            i += 1
    ax.bar(bars_x, bars_y, color=bars_c, edgecolor="black", linewidth=0.4)
    ax.set_xticks(bars_x)
    ax.set_xticklabels(labels)
    ax.set_ylabel("dG previsto (kcal/mol)")
    ax.axhline(0, color="black", lw=0.6)
    export(fig, "fig4_prodigy_affinity",
           {"method": "PRODIGY-prot 2.3.0, --selection A B, 25C, sobre a pose Boltz-2 model_0 protonada pH 8,2"})

# ---------- Fig 5: Interface overlap vs baseline MD (Boltz-2 only) ----------
baseline_top = {55, 56, 57, 58, 63, 99, 120, 121, 122, 123, 166, 167, 170, 175, 197, 218, 219, 220, 240, 242}
with style_context("nature", palette_name="okabe_ito_on_white") as info:
    w, h = mm_figsize(75, 65)
    fig, ax = plt.subplots(figsize=(w, h), layout="constrained")
    palette = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    bars_x, bars_y, bars_c, labels = [], [], [], []
    i = 0
    for sysname in SYS_LABEL:
        rows = list(csv.DictReader(open(BASE / sysname / "interacao" / "boltz_interface.csv")))
        resnums = {int(r["resnum"]) for r in rows}
        overlap = len(resnums & baseline_top)
        bars_x.append(i)
        bars_y.append(overlap)
        bars_c.append(palette[COLOR_KEY[sysname]])
        labels.append(SYS_LABEL[sysname].split(" ")[0])
        i += 1
    ax.bar(bars_x, bars_y, color=bars_c, edgecolor="black", linewidth=0.4)
    ax.axhline(len(baseline_top), color="0.4", ls="--", lw=0.7,
               label=f"total no baseline MD (n={len(baseline_top)})")
    ax.set_xticks(bars_x)
    ax.set_xticklabels(labels)
    ax.set_ylabel("residuos em comum c/ MD 100 ns")
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.32), frameon=False)
    export(fig, "fig5_interface_overlap",
           {"reference": "DN2954-GORE12T/MD/dn2954-gore12t/analise_extra/interface_residues.csv, top-20 por max_contact_freq"})

print("\nTodas as figuras exportadas em", OUT)
