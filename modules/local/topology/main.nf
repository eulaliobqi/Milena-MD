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
    // Força de campo externa (ex.: charmm36*, não empacotada no GROMACS stock):
    // pdb2gmx procura <nome>.ff primeiro no cwd -- link simbólico resolve sem
    // copiar os ~18 MB do .ff a cada work dir.
    def ff_link = params.forcefield_dir
        ? "ln -sfn ${params.forcefield_dir} ./${new File(params.forcefield_dir).name}"
        : ''
    """
    ${ff_link}
    printf '0\\n0\\n0\\n0\\n' | ${params.gmx_cmd} pdb2gmx \\
        -f ${complexo_pdb} \\
        -o complexo.gro -p topol.top -i posre.itp \\
        -ff ${params.forcefield} -water ${params.water} \\
        -ignh -ter -chainsep ter -merge no \\
        2>&1 | tee pdb2gmx.log
    """
}
