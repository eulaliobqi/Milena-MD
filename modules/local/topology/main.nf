process TOPOLOGY {
    tag "${meta.id}"
    label 'process_medium'

    publishDir { "${params.outdir}/${meta.id}/topo" }, mode: 'copy'

    input:
    tuple val(meta), path(complexo_pdb)

    output:
    tuple val(meta), path("complexo.gro"), path("topol.top"), path("*.itp"), emit: topology

    script:
    // 2 cadeias (A=receptor, B=ligante) → 4 seleções de terminais (N/C por cadeia)
    """
    printf '0\\n0\\n0\\n0\\n' | ${params.gmx_cmd} pdb2gmx \\
        -f ${complexo_pdb} \\
        -o complexo.gro -p topol.top -i posre.itp \\
        -ff ${params.forcefield} -water ${params.water} \\
        -ignh -ter -chainsep ter -merge no \\
        2>&1 | tee pdb2gmx.log
    """
}
