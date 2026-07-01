process STABILITY_FILTER {
    tag "${meta.id}"
    label 'process_low'
    errorStrategy 'ignore'

    publishDir { "${params.outdir}/${meta.id}/mmgbsa_prep" }, mode: 'copy'

    input:
    tuple val(meta), path(md_tpr), path(md_xtc), path(xvg_files), path(lig_ndx)

    output:
    tuple val(meta), path(md_tpr), path("stable.xtc"), path(lig_ndx), emit: stable
    tuple val(meta), path("stability_report.txt"),                     emit: report

    script:
    def min_frac  = params.stability_min_frac   ?: 0.3
    def window    = params.stability_window_ns  ?: 5.0
    def threshold = params.stability_threshold  ?: 0.15
    """
    echo "=== STABILITY FILTER: ${meta.id} ===" >&2

    # Encontra rmsd_backbone.xvg entre os arquivos encenados
    RMSD_FILE=\$(ls rmsd_backbone.xvg 2>/dev/null | head -1)

    if [ -n "\${RMSD_FILE}" ]; then
        START_PS=\$(stability_filter.py \\
            --rmsd-xvg "\${RMSD_FILE}" \\
            --window-ns   ${window} \\
            --sd-threshold ${threshold} \\
            --min-frac    ${min_frac} \\
            --report stability_report.txt)
        echo "[STABILITY] Plateau detectado — início estável: \${START_PS} ps" >&2
    else
        echo "AVISO: rmsd_backbone.xvg não encontrado — usando min_frac=${min_frac}" >&2
        START_PS=\$(python3 -c "print(int(${params.time_ns} * 1000 * ${min_frac}))")
        {
            echo "=== Relatório de Estabilidade ==="
            echo "Trajetória total        : ${params.time_ns}.0 ns"
            echo "Início fase estável     : \${START_PS} ps"
            echo "Método                  : fallback (rmsd_backbone.xvg ausente)"
        } > stability_report.txt
    fi

    # Extrai porção estável da trajetória
    echo "[STABILITY] Extraindo trajetória estável a partir de \${START_PS} ps" >&2
    echo "System" | ${params.gmx_cmd} trjconv \\
        -s ${md_tpr} \\
        -f ${md_xtc} \\
        -o stable.xtc \\
        -b \${START_PS} \\
        2>&1 | tail -5

    # Reporta frames resultantes
    N_FRAMES=\$(${params.gmx_cmd} check -f stable.xtc 2>&1 \\
        | grep -E "^Last frame" | awk '{print \$3}' || echo "?")
    echo "[OK] stable.xtc gerado com \${N_FRAMES} frames" | tee -a stability_report.txt
    """
}
