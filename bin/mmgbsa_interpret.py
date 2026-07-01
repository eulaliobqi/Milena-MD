#!/usr/bin/env python3
"""
Interpretação automática dos resultados MM-GBSA para complexo
proteína–peptídeo (tripsina digestiva + inibidor).

Saídas
------
mmgbsa_summary.txt   — valores-chave em texto legível
painel_mmgbsa.png    — painel 3×2 (ΔG, distribuição, componentes,
                        convergência, decomposição, hotspots)

Uso
---
mmgbsa_interpret.py --mmgbsa-csv mmgbsa_results.csv \
    [--decomp-csv decomp_results.csv] \
    [--stability-report stability_report.txt] \
    [--titulo "Sistema X"] \
    [--output-panel painel_mmgbsa.png] \
    [--output-summary mmgbsa_summary.txt]
"""
import argparse, os, sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.rcParams['font.size'] = 10
mpl.rcParams['axes.titlesize'] = 11


# ── Leitores ──────────────────────────────────────────────────────────────────

def load_mmgbsa_csv(path):
    """
    Lê CSV do gmx_MMPBSA. Retorna dict {colname: np.array} ou None.
    Tolerante a diferentes versões do gmx_MMPBSA (nomes de colunas variáveis).
    """
    if not path or not os.path.exists(path):
        return None
    try:
        with open(path) as f:
            lines = [l.rstrip() for l in f
                     if l.strip() and not l.startswith('#')]
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
        return {h: arr[:, i] for i, h in enumerate(header)
                if i < arr.shape[1]}
    except Exception as e:
        sys.stderr.write(f"[AVISO] Erro ao ler {path}: {e}\n")
        return None


def load_decomp_csv(path):
    """
    Lê CSV de decomposição por resíduo do gmx_MMPBSA.
    Retorna dict com chaves 'res_id', 'res_name', 'total' (arrays).
    """
    if not path or not os.path.exists(path):
        return None
    try:
        with open(path) as f:
            lines = [l.rstrip() for l in f
                     if l.strip() and not l.startswith('#')]
        if len(lines) < 2:
            return None
        header = [h.strip() for h in lines[0].split(',')]
        res_ids, res_names, rows = [], [], []
        for ln in lines[1:]:
            parts = [p.strip() for p in ln.split(',')]
            if len(parts) < 3:
                continue
            try:
                res_ids.append(int(parts[0]))
                res_names.append(parts[1] if len(parts) > 1 else '?')
                rows.append([float(x) for x in parts[2:] if x])
            except (ValueError, IndexError):
                pass
        if not rows:
            return None
        max_cols = max(len(r) for r in rows)
        arr = np.array([r + [np.nan] * (max_cols - len(r)) for r in rows])
        # Coluna TOTAL é a última coluna numérica
        total = np.nansum(arr, axis=1) if arr.shape[1] > 1 else arr[:, 0]
        return {
            'res_id':   np.array(res_ids),
            'res_name': res_names,
            'total':    total,
            'data':     arr,
            'header':   header[2:],
        }
    except Exception as e:
        sys.stderr.write(f"[AVISO] Erro ao ler decomp {path}: {e}\n")
        return None


# ── Interpretação biológica ────────────────────────────────────────────────────

def find_total_column(mmgbsa):
    """Retorna o nome da coluna de ΔG total."""
    candidates = ['TOTAL', 'DELTA total', 'delta_total', 'DeltaG',
                  'DELTA TOTAL', 'delta total', 'DELTE_TOTAL']
    for k in candidates:
        if k in mmgbsa:
            return k
    # Fallback: última coluna
    return list(mmgbsa.keys())[-1]


def classify_binding(dg_mean):
    """Classifica a qualidade de ligação para um peptídeo inibidor."""
    if dg_mean < -15:
        return "FORTE INIBIDOR",     "darkgreen"
    if dg_mean < -8:
        return "INIBIDOR MODERADO",  "forestgreen"
    if dg_mean < -3:
        return "INIBIDOR FRACO",     "darkorange"
    if dg_mean < 0:
        return "LIGAÇÃO MARGINAL",   "orange"
    return "SEM LIGAÇÃO DETECTADA",  "red"


# ── Plot helpers ──────────────────────────────────────────────────────────────

def _missing(ax, msg="(dados não encontrados)"):
    ax.text(0.5, 0.5, msg, ha='center', va='center',
            transform=ax.transAxes, fontsize=9, color='gray')
    ax.set_xticks([]); ax.set_yticks([])


def plot_dg_timeseries(ax, total, dg_mean, dg_std, classi, cor):
    frames = np.arange(len(total))
    ax.plot(frames, total, lw=0.7, color='darkred', alpha=0.65)
    ax.axhline(dg_mean, ls='--', color='red', lw=1.4,
               label=f'Média: {dg_mean:.2f} kcal/mol')
    ax.fill_between(frames,
                    dg_mean - dg_std, dg_mean + dg_std,
                    alpha=0.12, color='red', label=f'±1 SD={dg_std:.2f}')
    ax.set_xlabel("Frame")
    ax.set_ylabel("ΔG bind (kcal/mol)")
    ax.set_title(f"Energia de Ligação MM-GBSA\n{classi}")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)


def plot_dg_histogram(ax, total, dg_mean, cor):
    ax.hist(total, bins=min(30, len(total)//3 + 1),
            color=cor, alpha=0.75, edgecolor='black', lw=0.4)
    ax.axvline(dg_mean, color='black', lw=1.5, ls='--',
               label=f'μ = {dg_mean:.2f}')
    ax.set_xlabel("ΔG bind (kcal/mol)")
    ax.set_ylabel("Frequência")
    ax.set_title("Distribuição de ΔG bind")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)


def plot_components(ax, mmgbsa):
    mapping = [
        ('VDWAALS', 'VdW',         'steelblue'),
        ('EEL',     'Eletrost.',   'coral'),
        ('EGB',     'GB solvat.',  'mediumpurple'),
        ('ESURF',   'SA (np)',     'gold'),
    ]
    avail = [(l, c, mmgbsa[k]) for k, l, c in mapping if k in mmgbsa]
    if not avail:
        _missing(ax, "(componentes não encontrados)")
        return
    labels = [x[0] for x in avail]
    means  = [x[2].mean() for x in avail]
    stds   = [x[2].std()  for x in avail]
    colors = [x[1] for x in avail]
    bars = ax.bar(labels, means, yerr=stds, color=colors, alpha=0.85,
                  capsize=5, edgecolor='black', lw=0.5)
    ax.axhline(0, color='black', lw=0.8)
    ax.set_ylabel("Energia (kcal/mol)")
    ax.set_title("Componentes MM-GBSA (média ± DP)")
    for bar, m in zip(bars, means):
        offset = 2.0 if m >= 0 else -4.0
        ax.text(bar.get_x() + bar.get_width() / 2,
                m + offset, f"{m:.1f}",
                ha='center', fontsize=8, fontweight='bold')
    ax.grid(alpha=0.3, axis='y')


def plot_convergence(ax, total, dg_mean):
    cum_mean = np.cumsum(total) / (np.arange(len(total)) + 1)
    cum_std  = np.array([total[:i+1].std() for i in range(len(total))])
    frames   = np.arange(len(total))
    ax.plot(frames, cum_mean, color='darkblue', lw=1.2,
            label='Média cumulativa')
    ax.fill_between(frames,
                    cum_mean - cum_std, cum_mean + cum_std,
                    alpha=0.12, color='blue')
    ax.axhline(dg_mean, ls='--', color='red', lw=0.9, alpha=0.6,
               label=f'ΔG = {dg_mean:.2f}')
    ax.set_xlabel("Frame")
    ax.set_ylabel("ΔG bind (kcal/mol)")
    ax.set_title("Convergência de ΔG bind")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)


def plot_decomp_bar(ax, decomp):
    if decomp is None:
        _missing(ax, "(decomposição não disponível)")
        return
    totals = decomp['total']
    colors = ['darkred' if v < 0 else 'steelblue' for v in totals]
    ax.bar(decomp['res_id'], totals, color=colors, alpha=0.8, width=0.9)
    ax.axhline(0, color='black', lw=0.8)
    ax.set_xlabel("Resíduo")
    ax.set_ylabel("ΔG contribuição (kcal/mol)")
    ax.set_title("Decomposição por Resíduo\n"
                 "(vermelho = favorável, azul = desfavorável)")
    ax.grid(alpha=0.3, axis='y')


def plot_hotspots(ax, decomp, top_n=10):
    if decomp is None:
        _missing(ax, "(hotspots não disponíveis)")
        return
    totals = decomp['total']
    order  = np.argsort(totals)
    idx    = order[:min(top_n, len(order))]
    labels = [f"{decomp['res_name'][i]}{decomp['res_id'][i]}" for i in idx]
    vals   = totals[idx]
    colors = ['darkred' if v < 0 else 'navy' for v in vals]
    ax.barh(range(len(idx)), vals, color=colors, alpha=0.85)
    ax.set_yticks(range(len(idx)))
    ax.set_yticklabels(labels, fontsize=9)
    ax.axvline(0, color='black', lw=0.8)
    ax.set_xlabel("ΔG contribuição (kcal/mol)")
    ax.set_title(f"Top-{top_n} Hotspots Energéticos")
    ax.grid(alpha=0.3, axis='x')


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--mmgbsa-csv',       required=True)
    ap.add_argument('--decomp-csv',       default=None)
    ap.add_argument('--stability-report', default=None)
    ap.add_argument('--titulo',           default='MM-GBSA — Tripsina + Peptídeo')
    ap.add_argument('--output-panel',     default='painel_mmgbsa.png')
    ap.add_argument('--output-summary',   default='mmgbsa_summary.txt')
    args = ap.parse_args()

    mmgbsa = load_mmgbsa_csv(args.mmgbsa_csv)
    decomp = load_decomp_csv(args.decomp_csv) if args.decomp_csv else None

    if mmgbsa is None:
        sys.stderr.write("ERRO: não foi possível ler mmgbsa_results.csv\n")
        sys.exit(1)

    total_key = find_total_column(mmgbsa)
    total     = mmgbsa[total_key]
    dg_mean   = total.mean()
    dg_std    = total.std()
    dg_sem    = dg_std / max(1, np.sqrt(len(total)))
    classi, cor = classify_binding(dg_mean)

    # ── Texto resumo ──────────────────────────────────────────────────────────
    lines = [
        f"=== Resultado MM-GBSA: {args.titulo} ===",
        "",
        f"ΔGbind (média)   : {dg_mean:.3f} kcal/mol",
        f"ΔGbind (±SD)     : ±{dg_std:.3f} kcal/mol",
        f"ΔGbind (±SEM)    : ±{dg_sem:.3f} kcal/mol",
        f"N frames          : {len(total)}",
        f"Classificação     : {classi}",
        "",
        "Componentes energéticas (média ± DP):",
    ]
    comp_map = {
        'VDWAALS':        'Van der Waals',
        'EEL':            'Eletrostática',
        'EGB':            'Solvat. GB',
        'ESURF':          'SASA não-polar',
        'DELTA G gas':    'ΔG fase gasosa',
        'DELTA G solv':   'ΔG solvat.',
    }
    for k, label in comp_map.items():
        if k in mmgbsa:
            v = mmgbsa[k]
            lines.append(f"  {label:22s}: {v.mean():8.3f} ± {v.std():.3f} kcal/mol")

    lines += ["", "Interpretação biológica:"]
    if dg_mean < -5:
        lines.append("  Evidência de ligação significativa. Potencial bioinseticida.")
        lines.append("  Peptídeo LALAK mostra afinidade pela tripsina digestiva.")
    elif dg_mean < 0:
        lines.append("  Ligação fraca ou marginal. Otimização estrutural recomendada.")
    else:
        lines.append("  Sem evidência de ligação estável. Revisar pose e parâmetros.")

    if decomp is not None:
        hotspot_mask = np.abs(decomp['total']) > 0.5
        n_hot = hotspot_mask.sum()
        lines += ["", f"Hotspots (|ΔG| > 0.5 kcal/mol): {n_hot} resíduos"]
        order = np.argsort(decomp['total'])
        for i in order[:10]:
            if abs(decomp['total'][i]) > 0.5:
                lines.append(f"  {decomp['res_name'][i]}{decomp['res_id'][i]:4d}"
                              f"  {decomp['total'][i]:+.3f} kcal/mol")

    # Lê estabilidade se disponível
    if args.stability_report and os.path.exists(args.stability_report):
        lines += ["", "--- Estabilidade da trajetória ---"]
        with open(args.stability_report) as f:
            lines += [l.rstrip() for l in f]

    with open(args.output_summary, 'w') as fh:
        fh.write('\n'.join(lines) + '\n')

    for ln in lines:
        print(ln)

    # ── Painel gráfico ────────────────────────────────────────────────────────
    has_decomp = decomp is not None
    nrows = 3 if has_decomp else 2
    fig, axes = plt.subplots(nrows, 2, figsize=(14, nrows * 4.5))
    fig.suptitle(args.titulo, fontsize=13, fontweight='bold')

    plot_dg_timeseries(axes[0, 0], total, dg_mean, dg_std, classi, cor)
    plot_dg_histogram( axes[0, 1], total, dg_mean, cor)
    plot_components(   axes[1, 0], mmgbsa)
    plot_convergence(  axes[1, 1], total, dg_mean)

    if has_decomp:
        plot_decomp_bar(axes[2, 0], decomp)
        plot_hotspots(  axes[2, 1], decomp)

    plt.tight_layout()
    plt.savefig(args.output_panel, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"\nPainel salvo: {args.output_panel}")


if __name__ == '__main__':
    main()
