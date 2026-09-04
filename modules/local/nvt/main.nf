process NVT {
    tag "${meta.id}"
    label 'process_gpu'

    publishDir { "${params.outdir}/${meta.id}/nvt" }, mode: 'copy'

    input:
    tuple val(meta), path(em_gro), path(top, stageAs: 'input.top'), path(itps, stageAs: 'itp_in/*')

    output:
    tuple val(meta), path("nvt.gro"), path("nvt.cpt"), path("topol.top"), path("*.itp"), emit: system

    script:
    def gpu_flags = params.use_gpu ? "-nb gpu -pme gpu -bonded gpu -gpu_id ${params.gpu_id}" : ""
    def mpi       = params.mpi_cmd ?: ""
    def temp      = params.temperature
    // CHARMM36 exige vdW com switching (Force-switch); Cut-off simples é
    // o protocolo recomendado só p/ AMBER (ver minimization/main.nf).
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

    cat > nvt.mdp << MDP_EOF
define          = -DPOSRES
integrator      = md
dt              = 0.002
nsteps          = 100000
nstxout         = 1000
nstvout         = 1000
nstenergy       = 500
cutoff-scheme   = Verlet
nstlist         = 20
coulombtype     = PME
rcoulomb        = 1.2
${vdw_block}
constraints     = h-bonds
constraint-algorithm = LINCS
continuation    = no
gen-vel         = yes
gen-temp        = ${temp}
gen-seed        = -1
tcoupl          = V-rescale
tc-grps         = Protein Non-Protein
tau-t           = 0.1 0.1
ref-t           = ${temp} ${temp}
pcoupl          = no
pbc             = xyz
MDP_EOF

    ${params.gmx_cmd} grompp \\
        -f nvt.mdp \\
        -c ${em_gro} -r ${em_gro} \\
        -p topol.top -o nvt.tpr \\
        -maxwarn ${params.maxwarn}

    ${mpi} ${params.gmx_cmd} mdrun \\
        -v -deffnm nvt \\
        -ntomp ${params.ntomp} \\
        -pin on ${gpu_flags}
    """
}
