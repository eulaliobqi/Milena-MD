process PROLIF_FINGERPRINT {
    tag "${meta.id}"
    label 'process_medium'
    errorStrategy 'ignore'

    publishDir { "${params.outdir}/${meta.id}/analise_extra" }, mode: 'copy'

    input:
    tuple val(meta), path(complexo_pdb), path(md_tpr), path(stable_xtc)

    output:
    tuple val(meta), path("prolif_fingerprint.csv"), path("prolif_log.txt"),
                     path("prolif_heatmap.png", optional: true), emit: results

    script:
    """
    echo "=== PROLIF_FINGERPRINT: ${meta.id} ===" >&2

    prolif_fingerprint.py \\
        --complexo-pdb ${complexo_pdb} \\
        --tpr ${md_tpr} \\
        --xtc ${stable_xtc} \\
        --out-dir .

    echo "[OK] ProLIF (ou fallback) concluido para ${meta.id}" >&2
    """
}
