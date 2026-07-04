process PHARMACOPHORE_PROFILE {
    tag "${meta.id}"
    label 'process_low'
    errorStrategy 'ignore'

    publishDir { "${params.outdir}/${meta.id}/analise_extra" }, mode: 'copy'

    input:
    tuple val(meta), path(interface_residues_csv)

    output:
    tuple val(meta), path("pharmacophore_profile.csv"),
                     path("pharmacophore_summary.png"), emit: results

    script:
    """
    echo "=== PHARMACOPHORE_PROFILE: ${meta.id} ===" >&2

    pharmacophore_profile.py \\
        --interface-csv ${interface_residues_csv} \\
        --threshold 0.2 \\
        --out-dir .

    echo "[OK] Perfil farmacoforico concluido para ${meta.id}" >&2
    """
}
