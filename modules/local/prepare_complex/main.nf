process PREPARE_COMPLEX {
    tag "${meta.id}"
    label 'process_low'

    publishDir { "${params.outdir}/${meta.id}/prep" }, mode: 'copy'

    input:
    tuple val(meta), path(receptor), path(ligand)

    output:
    tuple val(meta), path("complexo.pdb"), emit: complexo

    script:
    // ff-style só importa se o receptor ainda tiver HIS/CYS "cru" (uso
    // standalone do script) -- no pipeline normal, PREPARE_PH já renomeou
    // tudo via pdb2pqr_process.py; passar aqui de qualquer forma mantém os
    // dois pontos consistentes e cobre o caso standalone do script.
    def ff_style = params.forcefield.toString().startsWith('charmm36') ? 'charmm' : 'amber'
    """
    prepare_complex.py \\
        --receptor ${receptor} \\
        --ligante  ${ligand} \\
        --pH ${params.pH} \\
        --ff-style ${ff_style} \\
        --out-dir .
    """
}
