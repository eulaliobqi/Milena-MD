process MINIMIZATION {
    tag "${meta.id}"
    label 'process_gpu'

    publishDir { "${params.outdir}/${meta.id}/em" }, mode: 'copy'

    input:
    tuple val(meta), path(ions_gro), path(top, stageAs: 'input.top'), path(itps, stageAs: 'itp_in/*')

    output:
    tuple val(meta), path("em.gro"), path("topol.top"), path("*.itp"), emit: system

    script:
    // steep NÃO suporta PME GPU → usa apenas -nb gpu (sem -pme gpu)
    def gpu_flags = params.use_gpu ? "-nb gpu -gpu_id ${params.gpu_id}" : ""
    def mpi       = params.mpi_cmd  ?: ""
    // CHARMM36 exige vdW com switching (Force-switch, rvdw-switch=1.0-1.2 nm,
    // sem DispCorr) -- Cut-off simples é o protocolo recomendado só p/ AMBER.
    def is_charmm = params.forcefield.toString().startsWith('charmm36')
    def vdw_block = is_charmm
        ? "vdwtype         = Cut-off\nvdw-modifier    = Force-switch\nrvdw-switch     = 1.0\nrvdw            = 1.2\nDispCorr        = no"
        : "vdwtype         = Cut-off\nrvdw            = 1.2"
    // topol.top inclui a .ff externa por caminho relativo (ver box_solvate_ions) --
    // precisa dela de novo aqui, senão grompp falha em resolver o #include.
    def ff_link = params.forcefield_dir
        ? "ln -sfn ${params.forcefield_dir} ./${new File(params.forcefield_dir).name}"
        : ''
    """
    ${ff_link}
    cp ${top} topol.top
    cp itp_in/*.itp .

    # Estágio 1: steepest descent até emtol frouxo ou nsteps -- remove clashes
    # graves rápido, mas converge devagar perto do mínimo.
    cat > em_steep.mdp << MDP_EOF
integrator      = steep
emtol           = 1000.0
emstep          = 0.01
nsteps          = 50000
cutoff-scheme   = Verlet
nstlist         = 10
coulombtype     = PME
rcoulomb        = 1.2
${vdw_block}
pbc             = xyz
MDP_EOF

    ${params.gmx_cmd} grompp \\
        -f em_steep.mdp -c ${ions_gro} \\
        -p topol.top -o em_steep.tpr \\
        -maxwarn ${params.maxwarn}

    ${mpi} ${params.gmx_cmd} mdrun \\
        -v -deffnm em_steep \\
        -ntomp ${params.ntomp} \\
        -pin on ${gpu_flags}

    # Estágio 2: segunda passada de steep a partir do estágio 1, emtol mais
    # apertado -- 'cg' foi tentado aqui e FALHOU em produção ("The
    # coordinates could not be constrained. Minimizer 'cg' can not handle
    # constraint failures") -- limitação conhecida do GROMACS: 'cg' não sabe
    # rejeitar/reduzir passo quando SETTLE (água) falha a restringir, ao
    # contrário de 'steep'; como o sistema é sempre solvatado (água sempre
    # com SETTLE), 'cg' nunca é seguro aqui. 'steep' de novo com emtol mais
    # baixo ainda entrega o refinamento adicional pretendido, sem essa
    # incompatibilidade.
    cat > em_cg.mdp << MDP_EOF
integrator      = steep
emtol           = 100.0
emstep          = 0.005
nsteps          = 50000
cutoff-scheme   = Verlet
nstlist         = 10
coulombtype     = PME
rcoulomb        = 1.2
${vdw_block}
pbc             = xyz
MDP_EOF

    ${params.gmx_cmd} grompp \\
        -f em_cg.mdp -c em_steep.gro \\
        -p topol.top -o em_cg.tpr \\
        -maxwarn ${params.maxwarn}

    ${mpi} ${params.gmx_cmd} mdrun \\
        -v -deffnm em_cg \\
        -ntomp ${params.ntomp} \\
        -pin on ${gpu_flags}

    cp em_cg.gro em.gro
    """
}
