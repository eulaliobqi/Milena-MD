#!/usr/bin/env python3
"""Fig. 6 -- pose Boltz-2 (model_0) adotada como estrutura inicial de MD,
para os dois sistemas do Bloco 1, renderizada em PyMOL 3.1.0 no servidor
(scripts/render_boltz_poses.py) e composta aqui com os helpers do skill
scientific-visualization (mesmo padrao de make_figures_pub.py).

Escopo corrigido em 2026-09-03: mostra so a pose que segue para a MD
(Boltz-2). A comparacao com a pose HADDOCK descartada (divergencia de
conformacao, investigada em comparativo_vs_DN2954-GORE12T.md secao 7)
fica arquivada em figuras/arquivo_divergencia_haddock_boltz/, fora da
numeracao principal do artigo.

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


img_a = autocrop(mpimg.imread(RENDER / "boltz_pose_gore12t.png"))
img_b = autocrop(mpimg.imread(RENDER / "boltz_pose_ggs3.png"))

ar = img_a.shape[1] / img_a.shape[0]  # mesma camera nos dois paineis
PANEL_W_MM = 85.0
PANEL_H_MM = PANEL_W_MM / ar
TITLE_MM = 6.0
LEGEND_MM = 12.0

fig_w_mm = 178.0
fig_h_mm = TITLE_MM + PANEL_H_MM + LEGEND_MM

with style_context("nature", palette_name="okabe_ito_on_white") as info:
    print("style:", info)
    fig, axes = plt.subplots(1, 2, figsize=mm_figsize(fig_w_mm, fig_h_mm), layout="constrained")

    for ax, img, panel, title in (
        (axes[0], img_a, "a", "DN2954 × GORE1-2T (21 aa)"),
        (axes[1], img_b, "b", "DN2954 × GORE1-2T(GGS)3 (75 aa)"),
    ):
        ax.imshow(img)
        ax.set_axis_off()
        ax.set_title(title, fontsize=7.5)
        ax.text(-0.02, 1.04, panel, transform=ax.transAxes, fontsize=8,
                fontweight="bold", va="bottom", ha="right")

    legend_elems = [
        Line2D([0], [0], color="#0072B2", lw=3, label="GORE1-2T (Boltz-2 model 0)"),
        Line2D([0], [0], color="#D55E00", lw=3, label="GORE1-2T(GGS)3 (Boltz-2 model 0)"),
        Line2D([0], [0], color="#009E73", lw=3, label="Tétrade catalítica (His79/Asp126/Ser222/Asp216)"),
        Line2D([0], [0], color="0.6", lw=3, label="Receptor DN2954 (tripsina)"),
    ]
    fig.legend(handles=legend_elems, loc="lower center", ncol=2, frameon=False,
               fontsize=7, bbox_to_anchor=(0.5, -0.02))

    report = export_figure(
        fig, OUT / "fig6_boltz_pose", formats=["pdf", "png"], dpi=600,
        bbox_inches=None, font_mode="truetype", overwrite=True,
        write_manifest=True,
        provenance={
            "raw_data": "novas-analises/{DN2954-GORE1-2T,DN2954-GORE1-2T-GGS3}/boltz2/, novas-analises/data/DN2954-receptor.pdb",
            "transformations": [
                "PyMOL 3.1.0 (scripts/render_boltz_poses.py, executado no ambiente structure do servidor): "
                "superposicao estrutural (cmd.align) da cadeia A (receptor) da pose Boltz-2 model_0 sobre "
                "data/DN2954-receptor.pdb; mesma orientacao/zoom/escala nos dois paineis",
                "composicao: autocrop de margem branca por bounding-box de conteudo (scripts/make_fig6_boltz_pose.py, "
                "funcao autocrop) -- recorte apenas, sem reamostragem/redimensionamento dos PNGs de origem",
                "nenhuma alteracao de coordenadas alem da superposicao rigida (sem minimizacao, sem edicao de pose)",
            ],
        },
    )
    plt.close(fig)
    for o in report["outputs"]:
        print(" ->", o["path"], o["size_bytes"], "bytes")
