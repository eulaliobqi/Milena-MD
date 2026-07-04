process FE_RERUN {
    tag "${meta.id}"
    label 'process_medium'
    errorStrategy 'ignore'

    publishDir { "${params.outdir}/${meta.id}/fe_estimate" }, mode: 'copy'

    input:
    tuple val(meta), path(complexo_gro), path(top, stageAs: 'input.top'),
                     path(itps, stageAs: 'itp_in/*'), path(mmgbsa_xtc), path(lig_ndx)

    output:
    tuple val(meta), path("interaction_energy.xvg"), emit: energy
    tuple val(meta), path("fe_rerun.log"),            emit: log

    script:
    """
    echo "=== FE_RERUN: ${meta.id} ===" >&2
    cp ${top} topol.top
    cp itp_in/*.itp .

    {
        echo "=== FE_RERUN: rerun com energygrps Receptor/Ligante ==="
        echo "Sistema: ${meta.id}"
        echo ""
    } > fe_rerun.log

    # rerun so recalcula energias sobre a trajetoria fornecida -- nsteps/dt
    # nao importam para -rerun, so os termos de energia/grupos precisam
    # bater com a producao original.
    cat > rerun.mdp << MDP_EOF
integrator           = md
dt                   = 0.002
nsteps               = 0
nstenergy            = 1
cutoff-scheme        = Verlet
nstlist              = 20
coulombtype          = PME
rcoulomb             = 1.2
vdwtype              = Cut-off
rvdw                 = 1.2
constraints          = h-bonds
constraint-algorithm = LINCS
pbc                  = xyz
energygrps           = Receptor Ligante
MDP_EOF

    # NOTA: Nextflow roda o script com "set -e" implicito -- cada comando
    # arriscado precisa de "|| true" para nao abortar o script inteiro antes
    # do fallback abaixo rodar (mesma licao do bug corrigido em MMGBSA_ROBUST:
    # exit != 0 com errorStrategy 'ignore' derruba TODOS os outputs, mesmo os
    # de fallback, se o script morrer antes de escreve-los).
    ${params.gmx_cmd} grompp \\
        -f rerun.mdp \\
        -c ${complexo_gro} \\
        -p topol.top \\
        -n ${lig_ndx} \\
        -o rerun.tpr \\
        -maxwarn ${params.maxwarn} \\
        2>&1 | tee -a fe_rerun.log || true

    if [ -f rerun.tpr ]; then
        ${params.mpi_cmd} ${params.gmx_cmd} mdrun \\
            -s rerun.tpr -rerun ${mmgbsa_xtc} \\
            -deffnm rerun -ntomp ${params.ntomp} \\
            2>&1 | tee -a fe_rerun.log || true
    else
        echo "ERRO: rerun.tpr nao gerado (grompp falhou)" | tee -a fe_rerun.log
    fi

    if [ -f rerun.edr ]; then
        printf 'Coul-SR:Receptor-Ligante\\nLJ-SR:Receptor-Ligante\\n0\\n' | \\
            ${params.gmx_cmd} energy -f rerun.edr -o interaction_energy.xvg \\
            2>&1 | tee -a fe_rerun.log || true
    fi

    if [ ! -f interaction_energy.xvg ]; then
        echo "ERRO: interaction_energy.xvg nao gerado" | tee -a fe_rerun.log
        echo "# rerun falhou -- ver fe_rerun.log" > interaction_energy.xvg
    else
        echo "[OK] interaction_energy.xvg gerado" | tee -a fe_rerun.log
    fi
    """
}
