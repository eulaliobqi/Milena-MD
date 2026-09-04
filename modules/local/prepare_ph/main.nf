process PREPARE_PH {
    tag "${meta.id}"
    label 'process_low'

    publishDir { "${params.outdir}/${meta.id}/prep_ph" }, mode: 'copy'

    input:
    tuple val(meta), path(receptor), path(ligand)

    output:
    tuple val(meta), path("receptor_ph.pdb"), path("ligand_ph.pdb"), emit: protonated
    tuple val(meta), path("*.propka"), optional: true, emit: propka

    script:
    def ph = params.pH
    // PDB2PQR só sabe emitir nomenclatura AMBER-flavored (--ff AMBER
    // --ffout AMBER, sempre, não tem --ffout CHARMM utilizável aqui) --
    // pdb2pqr_process.py é quem faz o remapeamento pra convenção da força
    // de campo final (HID/HIE/HIP+CYX pra amber, HSD/HSE/HSP+CYS2+ASPP+GLUP
    // pra charmm). Sem isso: grompp aborta com "Residue 'HIE' not found"
    // (confirmado em produção na 1a tentativa CHARMM36).
    def ff_style = params.forcefield.toString().startsWith('charmm36') ? 'charmm' : 'amber'
    """
    # ── Receptor ──────────────────────────────────────────────────────────────
    pdb2pqr --ff AMBER --ffout AMBER \\
        --titration-state-method propka --with-ph ${ph} \\
        --pdb-output receptor_raw.pdb \\
        --nodebump \\
        ${receptor} receptor.pqr

    pdb2pqr_process.py receptor_raw.pdb receptor_ph.pdb ${ff_style}

    # ── Ligante (peptídio) ─────────────────────────────────────────────────────
    # PDB2PQR pode falhar em peptídios curtos; usa o original como fallback
    if pdb2pqr --ff AMBER --ffout AMBER \\
        --titration-state-method propka --with-ph ${ph} \\
        --pdb-output ligand_raw.pdb \\
        --nodebump \\
        ${ligand} ligand.pqr ; then
        pdb2pqr_process.py ligand_raw.pdb ligand_ph.pdb ${ff_style}
    else
        echo "WARNING: pdb2pqr falhou no ligante, usando PDB original sem re-protonação"
        cp ${ligand} ligand_ph.pdb
    fi
    """
}
