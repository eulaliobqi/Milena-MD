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
    """
    cp ${top} topol.top
    cp itp_in/*.itp .

    cat > em.mdp << 'MDP_EOF'
integrator      = steep
emtol           = 1000.0
emstep          = 0.01
nsteps          = 50000
cutoff-scheme   = Verlet
nstlist         = 10
coulombtype     = PME
rcoulomb        = 1.2
vdwtype         = Cut-off
rvdw            = 1.2
pbc             = xyz
MDP_EOF

    ${params.gmx_cmd} grompp \\
        -f em.mdp -c ${ions_gro} \\
        -p topol.top -o em.tpr \\
        -maxwarn ${params.maxwarn}

    ${mpi} ${params.gmx_cmd} mdrun \\
        -v -deffnm em \\
        -ntomp ${params.ntomp} \\
        -pin on ${gpu_flags}
    """
}
