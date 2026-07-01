process MMGBSA_INTERPRET {
    tag "${meta.id}"
    label 'process_low'
    errorStrategy 'ignore'

    publishDir { "${params.outdir}/${meta.id}/mmgbsa" }, mode: 'copy'

    input:
    tuple val(meta),
        path(mmgbsa_csv),
        path(mmgbsa_dat),
        path(decomp_csv),
        path(stability_report)

    output:
    tuple val(meta),
        path("painel_mmgbsa.png"),
        path("mmgbsa_summary.txt"), emit: report

    script:
    def titulo = "${meta.id} — DM ${params.time_ns} ns @ pH ${params.pH} [MM-GBSA]"
    """
    # Ignora decomp vazia (apenas cabeçalho)
    DECOMP_ARG=""
    if [ -f ${decomp_csv} ] && [ \$(wc -l < ${decomp_csv}) -gt 1 ]; then
        DECOMP_ARG="--decomp-csv ${decomp_csv}"
    fi

    STAB_ARG=""
    [ -f ${stability_report} ] && STAB_ARG="--stability-report ${stability_report}"

    mmgbsa_interpret.py \\
        --mmgbsa-csv   ${mmgbsa_csv} \\
        \${DECOMP_ARG} \\
        \${STAB_ARG} \\
        --titulo       "${titulo}" \\
        --output-panel painel_mmgbsa.png \\
        --output-summary mmgbsa_summary.txt
    """
}
