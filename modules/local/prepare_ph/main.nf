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
    """
    # ── Receptor ──────────────────────────────────────────────────────────────
    pdb2pqr --ff AMBER --ffout AMBER \\
        --titration-state-method propka --with-ph ${ph} \\
        --pdb-output receptor_raw.pdb \\
        --nodebump \\
        ${receptor} receptor.pqr

    pdb2pqr_process.py receptor_raw.pdb receptor_ph.pdb

    # ── Ligante (peptídio) ─────────────────────────────────────────────────────
    # PDB2PQR pode falhar em peptídios curtos; usa o original como fallback
    if pdb2pqr --ff AMBER --ffout AMBER \\
        --titration-state-method propka --with-ph ${ph} \\
        --pdb-output ligand_raw.pdb \\
        --nodebump \\
        ${ligand} ligand.pqr ; then
        pdb2pqr_process.py ligand_raw.pdb ligand_ph.pdb
    else
        echo "WARNING: pdb2pqr falhou no ligante, usando PDB original sem re-protonação"
        cp ${ligand} ligand_ph.pdb
    fi
    """
}
