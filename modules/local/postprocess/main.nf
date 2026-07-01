process POSTPROCESS {
    tag "${meta.id}"
    label 'process_medium'

    publishDir { "${params.outdir}/${meta.id}/prod" }, mode: 'copy', overwrite: true

    input:
    tuple val(meta), path(md_tpr, stageAs: 'input.tpr'), path(md_xtc)

    output:
    tuple val(meta), path("md.tpr"), path("md_fit.xtc"), emit: fit

    script:
    """
    # 1) nojump: evita spikes de PBC na trajetória
    echo 0 | ${params.gmx_cmd} trjconv \\
        -s ${md_tpr} -f ${md_xtc} \\
        -o md_nojump.xtc -pbc nojump

    # 2) Corrige PBC: centra a proteína, compacta o solvente
    printf '1\\n0\\n' | ${params.gmx_cmd} trjconv \\
        -s ${md_tpr} -f md_nojump.xtc \\
        -o md_nopbc.xtc -pbc mol -center -ur compact

    # 3) Alinha pelo backbone (rot+trans fit)
    printf '4\\n0\\n' | ${params.gmx_cmd} trjconv \\
        -s ${md_tpr} -f md_nopbc.xtc \\
        -o md_fit.xtc -fit rot+trans

    # md.tpr é necessário nas análises — cria link local com nome fixo
    cp ${md_tpr} md.tpr

    rm -f md_nojump.xtc md_nopbc.xtc
    """
}
