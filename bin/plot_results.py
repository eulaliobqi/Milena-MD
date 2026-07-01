#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Gera painel de figuras a partir dos .xvg do GROMACS.

Todos os gráficos de série temporal mostram:
  - traço bruto (leve, transparente)
  - média móvel (linha sólida, janela configurável via --window-ns)
  - banda ±1 desvio padrão (área sombreada)

Uso:
  plot_results.py --analise-dir <dir> [--titulo T] [--window-ns N] [--mmgbsa-csv F]
"""
import argparse, os, sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.rcParams.update({
    'font.size':       10,
    'axes.titlesize':  11,
    'axes.labelsize':  10,
    'legend.fontsize':  8,
    'figure.dpi':     150,
})


# ── Leitores ──────────────────────────────────────────────────────────────────

def load_xvg(path):
    if not path or not os.path.exists(path):
        return None
    data = []
    for ln in open(path):
        if ln.startswith(('#', '@')):
            continue
        try:
            data.append([float(x) for x in ln.split()])
        except ValueError:
            pass
    return np.array(data) if data else None


def load_mmgbsa_csv(path):
    if not path or not os.path.exists(path):
        return None
    try:
        with open(path) as f:
            lines = [l.rstrip() for l in f if l.strip() and not l.startswith('#')]
        if len(lines) < 2:
            return None
        header = [h.strip() for h in lines[0].split(',')]
        rows = []
        for ln in lines[1:]:
            try:
                rows.append([float(x) for x in ln.split(',')])
            except ValueError:
                pass
        if not rows:
            return None
        arr = np.array(rows)
        return {h: arr[:, i] for i, h in enumerate(header) if i < arr.shape[1]}
    except Exception:
        return None


# ── Estatísticas de janela deslizante (numpy puro, O(n)) ──────────────────────

def rolling_stats(values, window):
    """Média móvel e desvio padrão por janela centrada. Sem dependências extras."""
    n   = len(values)
    w   = max(1, min(window, n))
    v   = values.astype(float)
    # padding com borda para não distorcer extremos
    pad  = w // 2
    vp   = np.pad(v, (pad, w - pad - 1), mode='edge')
    # cumsum para média eficiente
    cs   = np.cumsum(np.insert(vp, 0, 0.0))
    mean = (cs[w:] - cs[:-w]) / w
    cs2  = np.cumsum(np.insert(vp ** 2, 0, 0.0))
    var  = (cs2[w:] - cs2[:-w]) / w - mean ** 2
    std  = np.sqrt(np.maximum(var, 0.0))
    return mean[:n], std[:n]


# ── Helpers de plotagem ────────────────────────────────────────────────────────

def _missing(ax):
    ax.text(0.5, 0.5, "(file not found)", ha='center', va='center',
            transform=ax.transAxes, fontsize=9, color='gray')
    ax.set_xticks([]); ax.set_yticks([])


def plot_line(ax, data, ylabel, color, title=None, hline=None,
              xlabel="Time (ns)", window_ns=5.0):
    """Time series with rolling mean and ±1 SD band."""
    if data is None:
        _missing(ax); return

    t      = data[:, 0]
    values = data[:, 1]
    media  = values.mean()
    sdg    = values.std()

    dt     = float(t[1] - t[0]) if len(t) > 1 else 0.01
    window = max(3, int(window_ns / dt))

    rm, rs = rolling_stats(values, window)

    # traço bruto
    ax.plot(t, values, lw=0.4, color=color, alpha=0.20, zorder=1)
    # banda ±1 DP
    ax.fill_between(t, rm - rs, rm + rs, alpha=0.18, color=color, zorder=2)
    # média móvel
    ax.plot(t, rm, lw=1.6, color=color, zorder=3)

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    if title:
        ax.set_title(f"{title}\n(mean: {media:.3f} ± {sdg:.3f})")
    if hline is not None:
        ax.axhline(hline, ls='--', color='red', alpha=0.5, lw=0.9)
    ax.grid(alpha=0.25)


def plot_rmsf(ax, data):
    if data is None:
        _missing(ax); return
    ax.bar(data[:, 0], data[:, 1], width=1.0, color='seagreen', alpha=0.8)
    ax.bar(data[-5:, 0], data[-5:, 1], width=1.0, color='crimson',
           label='Ligand (last 5)')
    ax.legend()
    ax.set_xlabel("Residue")
    ax.set_ylabel("RMSF (nm)")
    ax.set_title("Per-residue fluctuation (RMSF)")
    ax.grid(alpha=0.25, axis='y')


def plot_sasa_dual(ax, sasa_prot, sasa_lig, window_ns):
    """SASA receptor (eixo esq.) + ligante (eixo dir.) com médias móveis."""
    has_p = sasa_prot is not None
    has_l = sasa_lig  is not None

    if not has_p and not has_l:
        _missing(ax); return

    cor_p, cor_l = 'steelblue', 'tomato'
    ax2 = ax.twinx()

    def _add(axis, data, color, label, ylabel):
        t, v  = data[:, 0], data[:, 1]
        dt    = float(t[1] - t[0]) if len(t) > 1 else 0.01
        win   = max(3, int(window_ns / dt))
        rm, rs = rolling_stats(v, win)
        axis.plot(t, v, lw=0.4, color=color, alpha=0.20)
        axis.fill_between(t, rm - rs, rm + rs, alpha=0.15, color=color)
        axis.plot(t, rm, lw=1.6, color=color, label=label)
        axis.set_ylabel(ylabel, color=color)
        axis.tick_params(axis='y', labelcolor=color)

    if has_p:
        _add(ax,  sasa_prot, cor_p, 'Receptor', 'SASA receptor (nm²)')
    if has_l:
        _add(ax2, sasa_lig,  cor_l, 'Ligand',   'SASA ligand (nm²)')

    ax.set_xlabel("Time (ns)")
    ax.set_title("Solvent Accessible Surface Area (SASA)")
    if has_p: ax.legend(loc='upper left')
    if has_l: ax2.legend(loc='upper right')
    ax.grid(alpha=0.25)


def plot_triad_lines(ax, xvg_dict, labels, window_ns):
    """4 distâncias à tríade/S1 com médias móveis + banda."""
    cores = ['forestgreen', 'royalblue', 'crimson', 'darkorange']
    keys  = ['dist_r1', 'dist_r2', 'dist_r3', 'dist_r4']
    any_  = False

    for key, label, cor in zip(keys, labels, cores):
        d = xvg_dict.get(key)
        if d is None:
            continue
        any_ = True
        t, v  = d[:, 0], d[:, 1]
        dt    = float(t[1] - t[0]) if len(t) > 1 else 0.01
        win   = max(3, int(window_ns / dt))
        rm, rs = rolling_stats(v, win)
        ax.plot(t, v, lw=0.3, color=cor, alpha=0.18)
        ax.fill_between(t, rm - rs, rm + rs, alpha=0.13, color=cor)
        ax.plot(t, rm, lw=1.5, color=cor, label=label)

    if not any_:
        _missing(ax); return

    ax.axhline(0.5, ls='--', color='gray', alpha=0.6, lw=0.9, label='0.5 nm')
    ax.set_xlabel("Time (ns)")
    ax.set_ylabel("Minimum distance (nm)")
    ax.set_title("Distances to Catalytic Residues")
    ax.legend()
    ax.grid(alpha=0.25)


def plot_triad_bars(ax, xvg_dict, labels):
    """Barras de média ± DP para os 4 resíduos."""
    cores = ['forestgreen', 'royalblue', 'crimson', 'darkorange']
    keys  = ['dist_r1', 'dist_r2', 'dist_r3', 'dist_r4']
    lbs, ms, ss, cs = [], [], [], []

    for key, label, cor in zip(keys, labels, cores):
        d = xvg_dict.get(key)
        if d is None:
            continue
        lbs.append(label)
        ms.append(d[:, 1].mean())
        ss.append(d[:, 1].std())
        cs.append(cor)

    if not lbs:
        _missing(ax); return

    bars = ax.bar(lbs, ms, yerr=ss, color=cs, alpha=0.80,
                  capsize=5, edgecolor='black', linewidth=0.5)
    ax.axhline(0.5, ls='--', color='gray', alpha=0.6, lw=0.9, label='0.5 nm')
    for bar, m in zip(bars, ms):
        ax.text(bar.get_x() + bar.get_width() / 2,
                bar.get_height() + 0.02,
                f"{m:.3f} nm", ha='center', fontsize=9)
    ax.set_ylabel("Mean distance (nm)")
    ax.set_title("Mean Catalytic Triad Occupancy\n(mean ± SD)")
    ax.legend()
    ax.grid(alpha=0.25, axis='y')


def plot_mmgbsa_total(ax, mmgbsa):
    total = None
    for key in ('TOTAL', 'DELTA total', 'delta_total', 'DeltaG'):
        if key in mmgbsa:
            total = mmgbsa[key]; break
    if total is None:
        _missing(ax); return
    frames = np.arange(len(total))
    media, std = total.mean(), total.std()
    ax.plot(frames, total, lw=0.5, color='darkred', alpha=0.35)
    rm, rs = rolling_stats(total, max(3, len(total) // 20))
    ax.fill_between(frames, rm - rs, rm + rs, alpha=0.15, color='darkred')
    ax.plot(frames, rm, lw=1.6, color='darkred')
    ax.axhline(media, ls='--', color='red', lw=1.2,
               label=f'Mean: {media:.2f} kcal/mol')
    ax.set_xlabel("Frame")
    ax.set_ylabel("ΔG bind (kcal/mol)")
    ax.set_title(f"MM-GBSA Binding Energy\n(mean: {media:.2f} ± {std:.2f} kcal/mol)")
    ax.legend(); ax.grid(alpha=0.25)


def plot_mmgbsa_components(ax, mmgbsa):
    components = {}
    for key, label in [('VDWAALS', 'VdW'), ('EEL', 'Elec'),
                       ('EGB', 'GB solv.'), ('ESURF', 'SA')]:
        if key in mmgbsa:
            components[label] = mmgbsa[key]
    if not components:
        _missing(ax); return
    labels = list(components.keys())
    means  = [v.mean() for v in components.values()]
    stds   = [v.std()  for v in components.values()]
    colors = ['steelblue', 'coral', 'mediumpurple', 'gold'][:len(labels)]
    bars = ax.bar(labels, means, yerr=stds, color=colors, alpha=0.85,
                  capsize=4, edgecolor='black', linewidth=0.5)
    ax.axhline(0, color='black', lw=0.8)
    ax.set_ylabel("Energy (kcal/mol)")
    ax.set_title("MM-GBSA Components\n(mean ± SD)")
    for bar, m in zip(bars, means):
        ax.text(bar.get_x() + bar.get_width() / 2,
                bar.get_height() + (2 if m >= 0 else -6),
                f"{m:.1f}", ha='center', fontsize=8)
    ax.grid(alpha=0.25, axis='y')


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--analise-dir", required=True)
    ap.add_argument("--titulo",    default="Molecular Dynamics")
    ap.add_argument("--window-ns", type=float, default=5.0,
                    help="Rolling window in ns (default: 5.0)")
    ap.add_argument("--mmgbsa-csv", default=None)
    ap.add_argument("--output",    default="painel_completo.png")
    args = ap.parse_args()

    D   = args.analise_dir
    wns = args.window_ns

    xvg = {
        "rmsd_bb":   load_xvg(os.path.join(D, "rmsd_backbone.xvg")),
        "rmsd_lig":  load_xvg(os.path.join(D, "rmsd_ligante.xvg")),
        "rmsf":      load_xvg(os.path.join(D, "rmsf_residuos.xvg")),
        "rg":        load_xvg(os.path.join(D, "gyrate.xvg")),
        "ncont":     load_xvg(os.path.join(D, "numcont.xvg")),
        "hbond":     load_xvg(os.path.join(D, "hbond.xvg")),
        "sasa_prot": load_xvg(os.path.join(D, "sasa_protein.xvg")),
        "sasa_lig":  load_xvg(os.path.join(D, "sasa_ligante.xvg")),
        "dist_r1":   load_xvg(os.path.join(D, "dist_r1.xvg")),
        "dist_r2":   load_xvg(os.path.join(D, "dist_r2.xvg")),
        "dist_r3":   load_xvg(os.path.join(D, "dist_r3.xvg")),
        "dist_r4":   load_xvg(os.path.join(D, "dist_r4.xvg")),
    }

    # Rótulos da tríade
    triad_info = os.path.join(D, "triad_info.txt")
    if os.path.exists(triad_info):
        raw = [l.strip() for l in open(triad_info) if l.strip()]
        triad_labels = (raw + ["?"] * 4)[:4]
    else:
        triad_labels = ["Res1", "Res2", "Res3", "S1"]

    mmgbsa     = load_mmgbsa_csv(args.mmgbsa_csv) if args.mmgbsa_csv else None
    has_mmgbsa = mmgbsa is not None
    has_sasa   = xvg["sasa_prot"] is not None or xvg["sasa_lig"] is not None
    has_triad  = any(xvg[k] is not None for k in ("dist_r1", "dist_r2",
                                                    "dist_r3", "dist_r4"))

    nrows = 3 + has_sasa + has_triad + (2 if has_mmgbsa else 0)
    fig, axes = plt.subplots(nrows, 2, figsize=(14, nrows * 4))
    fig.suptitle(args.titulo, fontsize=14, fontweight='bold')

    # ── Linha 1: RMSD ─────────────────────────────────────────────────────────
    media_bb = xvg["rmsd_bb"][:, 1].mean() if xvg["rmsd_bb"] is not None else None
    plot_line(axes[0, 0], xvg["rmsd_bb"],  "RMSD (nm)", "navy",
              "Backbone RMSD", hline=media_bb, window_ns=wns)
    plot_line(axes[0, 1], xvg["rmsd_lig"], "RMSD (nm)", "darkorange",
              "Ligand RMSD (peptide)", window_ns=wns)

    # ── Linha 2: RMSF + Rg ────────────────────────────────────────────────────
    plot_rmsf(axes[1, 0], xvg["rmsf"])
    plot_line(axes[1, 1], xvg["rg"],       "Rg (nm)",   "purple",
              "Radius of Gyration", window_ns=wns)

    # ── Linha 3: Contatos + H-bonds ───────────────────────────────────────────
    plot_line(axes[2, 0], xvg["ncont"],    "N. atoms", "teal",
              "Receptor–ligand contacts (<0.4 nm)", window_ns=wns)
    plot_line(axes[2, 1], xvg["hbond"],    "N. H-bonds", "indianred",
              "Receptor–ligand hydrogen bonds", window_ns=wns)

    next_row = 3

    # ── Linha SASA ────────────────────────────────────────────────────────────
    if has_sasa:
        plot_sasa_dual(axes[next_row, 0], xvg["sasa_prot"],
                       xvg["sasa_lig"], wns)
        plot_line(axes[next_row, 1], xvg["sasa_lig"], "SASA (nm²)", "tomato",
                  "Ligand SASA (burial)", window_ns=wns)
        next_row += 1

    # ── Linha Tríade ──────────────────────────────────────────────────────────
    if has_triad:
        plot_triad_lines(axes[next_row, 0], xvg, triad_labels, wns)
        plot_triad_bars( axes[next_row, 1], xvg, triad_labels)
        next_row += 1

    # ── Linha MM-GBSA ─────────────────────────────────────────────────────────
    if has_mmgbsa:
        plot_mmgbsa_total(     axes[next_row, 0], mmgbsa)
        plot_mmgbsa_components(axes[next_row, 1], mmgbsa)

    plt.tight_layout()
    out = os.path.join(D, args.output)
    plt.savefig(out, dpi=150, bbox_inches='tight')
    print(f"Saved: {out}")
    plt.close()

    # ── PNGs individuais ──────────────────────────────────────────────────────
    specs_line = [
        ("rmsd_bb",   "RMSD (nm)",    "navy",       "rmsd_bb.png"),
        ("rmsd_lig",  "RMSD (nm)",    "darkorange",  "rmsd_lig.png"),
        ("rg",        "Rg (nm)",      "purple",      "rg.png"),
        ("ncont",     "N. atoms",     "teal",        "ncont.png"),
        ("hbond",     "N. H-bonds",   "indianred",   "hbond.png"),
        ("sasa_lig",  "SASA (nm²)",   "tomato",      "sasa_ligante.png"),
        ("sasa_prot", "SASA (nm²)",   "steelblue",   "sasa_protein.png"),
    ]
    for key, ylabel, color, fname in specs_line:
        if xvg.get(key) is None:
            continue
        fig, ax = plt.subplots(figsize=(9, 5))
        plot_line(ax, xvg[key], ylabel, color, title=fname.replace('.png',''),
                  window_ns=wns)
        plt.tight_layout()
        plt.savefig(os.path.join(D, fname), dpi=150, bbox_inches='tight')
        plt.close()

    if xvg["rmsf"] is not None:
        fig, ax = plt.subplots(figsize=(9, 5))
        plot_rmsf(ax, xvg["rmsf"])
        plt.tight_layout()
        plt.savefig(os.path.join(D, "rmsf.png"), dpi=150, bbox_inches='tight')
        plt.close()

    if has_triad:
        fig, ax = plt.subplots(figsize=(9, 5))
        plot_triad_lines(ax, xvg, triad_labels, wns)
        plt.tight_layout()
        plt.savefig(os.path.join(D, "triad_distances.png"), dpi=150, bbox_inches='tight')
        plt.close()

    print("Done.")


if __name__ == "__main__":
    main()
