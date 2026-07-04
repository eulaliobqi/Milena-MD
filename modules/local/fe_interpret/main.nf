process FE_INTERPRET {
    tag "${meta.id}"
    label 'process_low'
    errorStrategy 'ignore'

    publishDir { "${params.outdir}/${meta.id}/fe_estimate" }, mode: 'copy'

    input:
    tuple val(meta), path(interaction_xvg)

    output:
    tuple val(meta), path("free_energy_estimate.txt"),
                     path("interaction_energy.png", optional: true), emit: results

    script:
    def titulo = "${meta.id} — Interaction Energy (Coul-SR + LJ-SR)"
    """
    echo "=== FE_INTERPRET: ${meta.id} ===" >&2

    interaction_entropy.py \\
        --xvg ${interaction_xvg} \\
        --temperature ${params.temperature} \\
        --titulo "${titulo}" \\
        --out-dir .

    echo "[OK] Estimativa de energia livre concluida para ${meta.id}" >&2
    """
}
