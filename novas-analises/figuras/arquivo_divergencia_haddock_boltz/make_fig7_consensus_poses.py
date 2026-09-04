#!/usr/bin/env python3
"""Fig. 7 -- poses consenso (HADDOCK3 padrao-ouro vs Boltz-2 fronteira,
model_0) para os dois sistemas do Bloco 1, renderizadas em PyMOL 3.1.0 no
servidor (scripts/render_consensus_poses.py) e compostas aqui com os
helpers do skill scientific-visualization (mesmo padrao de
make_figures_pub.py).

Paineis (a) e (b): regiao de interface (receptor inteiro + peptideo ate
15 A da superficie do receptor), mesma camera nos dois paineis (mesma
escala Angstrom/pixel -- os dois PNGs de origem so sao recortados aqui
para remover margem branca, nunca redimensionados).
Painel (c): pose HADDOCK completa do GORE1-2T(GGS)3 sem esse corte de
camera -- documenta que o melhor modelo HADDOCK para esse peptideo de 75 aa
manteve o conformero helice-alfa do ensemble de entrada quase reto ao
longo de toda a cadeia (~110 A), com apenas a extremidade N-terminal em
contato com o receptor; achado relevante para a preparacao da topologia de
MD (a cauda distal nao tem suporte experimental de conformacao e nao deve
ser tratada como pose validada).

Roda localmente (Windows), lendo os PNGs ja copiados de volta do servidor
em novas-analises/figuras/_render/.
"""
import sys
from pathlib import Path

import numpy as np

SKILL_SCRIPTS = Path(r"C:\Users\eulal\.claude\skills\scientific-visualization\scripts")
sys.path.insert(0, str(SKILL_SCRIPTS))

from style_presets import style_context  # noqa: E402
from figure_export import export_figure  # noqa: E402

import matplotlib.pyplot as plt  # noqa: E402
import matplotlib.image as mpimg  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

BASE = Path(__file__).resolve().parents[1]
RENDER = BASE / "figuras" / "_render"
OUT = BASE / "figuras"
OUT.mkdir(exist_ok=True)


def mm_figsize(width_mm, height_mm):
    return (width_mm / 25.4, height_mm / 25.4)


def autocrop(img, pad=25):
    """Recorta margem branca (opaca) ao redor do conteudo. So corta, nunca
    redimensiona -- preserva a escala Angstrom/pixel entre paineis que
    compartilham a mesma camera PyMOL."""
    rgb = img[..., :3]
    non_white = np.any(rgb < 0.98, axis=-1)
    rows = np.where(non_white.any(axis=1))[0]
    cols = np.where(non_white.any(axis=0))[0]
    r0, r1 = max(rows[0] - pad, 0), min(rows[-1] + pad, img.shape[0] - 1)
    c0, c1 = max(cols[0] - pad, 0), min(cols[-1] + pad, img.shape[1] - 1)
    return img[r0:r1 + 1, c0:c1 + 1]


img_a = autocrop(mpimg.imread(RENDER / "consensus_gore12t.png"))
img_b = autocrop(mpimg.imread(RENDER / "consensus_ggs3.png"))
img_c = autocrop(mpimg.imread(RENDER / "consensus_ggs3_full.png"))

ar_ab = img_a.shape[1] / img_a.shape[0]  # largura/altura (a e b usam a mesma camera)
ar_c = img_c.shape[1] / img_c.shape[0]

TOP_W_MM = 85.0   # largura de cada painel a/b
TOP_H_MM = TOP_W_MM / ar_ab
BOT_W_MM = 178.0  # largura do painel c
BOT_H_MM = BOT_W_MM / ar_c
TITLE_MM = 6.0
LEGEND_MM = 14.0
GAP_MM = 6.0

fig_w_mm = 178.0
fig_h_mm = TITLE_MM + TOP_H_MM + GAP_MM + TITLE_MM + BOT_H_MM + LEGEND_MM

with style_context("nature", palette_name="okabe_ito_on_white") as info:
    print("style:", info)
    fig = plt.figure(figsize=mm_figsize(fig_w_mm, fig_h_mm), layout="constrained")
    gs = fig.add_gridspec(
        2, 2,
        height_ratios=[TOP_H_MM + TITLE_MM, BOT_H_MM + TITLE_MM],
        hspace=0.12, wspace=0.04,
    )
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[1, :])

    for ax, img, panel, title in (
        (ax_a, img_a, "a", "DN2954 \u00d7 GORE1-2T (21 aa)"),
        (ax_b, img_b, "b", "DN2954 \u00d7 GORE1-2T(GGS)3 (75 aa)"),
    ):
        ax.imshow(img)
        ax.set_axis_off()
        ax.set_title(title, fontsize=7.5)
        ax.text(-0.02, 1.04, panel, transform=ax.transAxes, fontsize=8,
                fontweight="bold", va="bottom", ha="right")

    ax_c.imshow(img_c)
    ax_c.set_axis_off()
    ax_c.set_title(
        "Pose HADDOCK completa do GORE1-2T(GGS)3, mesma cena sem corte de camera",
        fontsize=7.5,
    )
    ax_c.text(-0.005, 1.04, "c", transform=ax_c.transAxes, fontsize=8,
              fontweight="bold", va="bottom", ha="right")

    legend_elems = [
        Line2D([0], [0], color="#0072B2", lw=3, label="Peptideo \u2014 pose HADDOCK3 (padr\u00e3o-ouro)"),
        Line2D([0], [0], color="#D55E00", lw=3, label="Peptideo \u2014 pose Boltz-2 model\u00a00 (fronteira)"),
        Line2D([0], [0], color="#009E73", lw=3, label="T\u00e9trade catal\u00edtica (His79/Asp126/Ser222/Asp216)"),
        Line2D([0], [0], color="0.6", lw=3, label="Receptor DN2954 (tripsina)"),
    ]
    fig.legend(handles=legend_elems, loc="lower center", ncol=2, frameon=False,
               fontsize=7, bbox_to_anchor=(0.5, -0.01))

    report = export_figure(
        fig, OUT / "fig7_consensus_poses", formats=["pdf", "png"], dpi=600,
        bbox_inches=None, font_mode="truetype", overwrite=True,
        write_manifest=True,
        provenance={
            "raw_data": "novas-analises/{DN2954-GORE1-2T,DN2954-GORE1-2T-GGS3}/{haddock,boltz2}/, novas-analises/data/DN2954-receptor.pdb",
            "transformations": [
                "PyMOL 3.1.0 (scripts/render_consensus_poses.py, executado no ambiente structure do servidor): "
                "superposicao estrutural (cmd.align) da cadeia A (receptor) de cada pose sobre data/DN2954-receptor.pdb; "
                "paineis a/b limitam a camera ao receptor inteiro + peptideo ate 15 A da superficie do receptor (mesma "
                "orientacao/zoom nos dois paineis, mesma escala Angstrom/pixel); painel c usa a mesma cena sem esse "
                "limite de camera",
                "composicao: autocrop de margem branca por bounding-box de conteudo (scripts/make_fig7_consensus_poses.py, "
                "funcao autocrop) -- recorte apenas, sem reamostragem/redimensionamento dos PNGs de origem",
                "nenhuma alteracao de coordenadas alem da superposicao rigida (sem minimizacao, sem edicao de pose)",
            ],
        },
    )
    plt.close(fig)
    for o in report["outputs"]:
        print(" ->", o["path"], o["size_bytes"], "bytes")
