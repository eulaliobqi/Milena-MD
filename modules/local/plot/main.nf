process PLOT {
    tag "${meta.id}"
    label 'process_low'

    publishDir { "${params.outdir}/${meta.id}/analise" }, mode: 'copy'

    input:
    tuple val(meta), path(xvg_files), path(lig_ndx), path(mmgbsa_csv)

    output:
    tuple val(meta), path("*.png"), emit: figures

    script:
    def titulo = "${meta.id} - DM ${params.time_ns} ns @ pH ${params.pH} (KCl ${params.nacl_conc}M)"
    """
    plot_results.py \\
        --analise-dir . \\
        --titulo "${titulo}" \\
        --window-ns ${params.window_ns} \\
        --mmgbsa-csv ${mmgbsa_csv}
    """
}
