process PREPARE_COMPLEX {
    tag "${meta.id}"
    label 'process_low'

    publishDir { "${params.outdir}/${meta.id}/prep" }, mode: 'copy'

    input:
    tuple val(meta), path(receptor), path(ligand)

    output:
    tuple val(meta), path("complexo.pdb"), emit: complexo

    script:
    """
    prepare_complex.py \\
        --receptor ${receptor} \\
        --ligante  ${ligand} \\
        --pH ${params.pH} \\
        --out-dir .
    """
}
