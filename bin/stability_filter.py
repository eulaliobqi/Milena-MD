#!/usr/bin/env python3
"""
Detecta o início da fase estável da trajetória a partir do RMSD do backbone.

Método: janela deslizante de desvio padrão. A fase estável começa quando o
SD local cai abaixo de (threshold_frac × SD_global) por pelo menos 3 janelas
consecutivas. Se a detecção falhar, usa os últimos (1 - min_frac) da trajetória.

Saída stdout : tempo de início em ps (inteiro), para capturar com $(...)
Arquivo      : stability_report.txt (criado no diretório corrente)
"""
import argparse, sys
import numpy as np


def load_xvg(path):
    data = []
    for ln in open(path):
        if ln.startswith(('#', '@')):
            continue
        try:
            vals = [float(x) for x in ln.split()]
            if len(vals) >= 2:
                data.append(vals[:2])
        except ValueError:
            pass
    return np.array(data) if data else None


def find_stable_start(times_ns, rmsd, window_ns=5.0, sd_threshold_frac=0.15,
                      min_frac=0.3, required_consecutive=3):
    """
    Detecta plateau do RMSD via janela deslizante de SD.

    Retorna o tempo (ns) a partir do qual a trajetória é considerada estável.
    Nunca retorna um valor anterior a (min_frac × duração_total).
    """
    n = len(times_ns)
    if n < 20:
        return times_ns[int(n * min_frac)]

    dt = (times_ns[-1] - times_ns[0]) / (n - 1)
    window = max(5, int(window_ns / dt))
    global_sd = rmsd.std()

    if global_sd < 1e-6:
        return times_ns[int(n * min_frac)]

    threshold = sd_threshold_frac * global_sd
    consecutive = 0
    plateau_start_idx = None

    for i in range(n - window):
        local_sd = rmsd[i:i + window].std()
        if local_sd < threshold:
            consecutive += 1
            if consecutive >= required_consecutive and plateau_start_idx is None:
                # Índice onde começa o bloco estável
                plateau_start_idx = max(0, i - (required_consecutive - 1) * window)
        else:
            consecutive = 0

    if plateau_start_idx is not None:
        start_ns = times_ns[plateau_start_idx]
    else:
        # Fallback: 40% inicial = pré-equilíbrio descartado
        start_ns = times_ns[int(n * min_frac)]

    # Garante que usamos pelo menos (1 - min_frac) da trajetória
    early_limit = times_ns[int(n * min_frac)]
    return max(start_ns, early_limit)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--rmsd-xvg',      required=True, help='Arquivo rmsd_backbone.xvg')
    ap.add_argument('--window-ns',     type=float, default=5.0,
                    help='Tamanho da janela deslizante em ns (default: 5.0)')
    ap.add_argument('--sd-threshold',  type=float, default=0.15,
                    help='Fração do SD global como limiar (default: 0.15)')
    ap.add_argument('--min-frac',      type=float, default=0.3,
                    help='Fração mínima descartada como equilíbrio (default: 0.3)')
    ap.add_argument('--report',        default='stability_report.txt')
    args = ap.parse_args()

    data = load_xvg(args.rmsd_xvg)
    if data is None or len(data) < 10:
        sys.stderr.write("AVISO: RMSD insuficiente — usando fallback 30%\n")
        print("0", flush=True)
        return

    times = data[:, 0]   # unidade: ns (GROMACS escreve -tu ns)
    rmsd  = data[:, 1]   # nm
    total_ns = times[-1]

    stable_start_ns = find_stable_start(
        times, rmsd,
        window_ns=args.window_ns,
        sd_threshold_frac=args.sd_threshold,
        min_frac=args.min_frac
    )
    stable_start_ps = int(round(stable_start_ns * 1000))

    stable_mask = times >= stable_start_ns
    rmsd_stable = rmsd[stable_mask]
    rmsd_all    = rmsd

    lines = [
        "=== Relatório de Estabilidade ===",
        f"Trajetória total        : {total_ns:.1f} ns",
        f"Início fase estável     : {stable_start_ns:.2f} ns  ({stable_start_ps} ps)",
        f"Duração fase estável    : {total_ns - stable_start_ns:.1f} ns  "
        f"({(total_ns - stable_start_ns) / total_ns * 100:.0f}% da trajetória)",
        f"RMSD médio (estável)    : {rmsd_stable.mean():.4f} ± {rmsd_stable.std():.4f} nm",
        f"RMSD médio (global)     : {rmsd_all.mean():.4f} ± {rmsd_all.std():.4f} nm",
        f"Frames estáveis         : {stable_mask.sum()} de {len(times)}",
    ]

    with open(args.report, 'w') as fh:
        fh.write('\n'.join(lines) + '\n')

    for ln in lines:
        sys.stderr.write(ln + '\n')

    # Apenas o valor numérico em stdout para captura pelo shell
    print(stable_start_ps, flush=True)


if __name__ == '__main__':
    main()
