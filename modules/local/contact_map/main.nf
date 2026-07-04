process CONTACT_MAP {
    tag "${meta.id}"
    label 'process_medium'
    errorStrategy 'ignore'

    publishDir { "${params.outdir}/${meta.id}/analise_extra" }, mode: 'copy'

    input:
    tuple val(meta), path(complexo_pdb), path(md_tpr), path(stable_xtc)

    output:
    tuple val(meta), path("contact_map.csv"), path("contact_map.png"),
                     path("interface_residues.csv"), emit: results

    script:
    """
    echo "=== CONTACT_MAP: ${meta.id} ===" >&2

    contact_map.py \\
        --complexo-pdb ${complexo_pdb} \\
        --tpr ${md_tpr} \\
        --xtc ${stable_xtc} \\
        --out-dir .

    echo "[OK] Mapa de contatos concluido para ${meta.id}" >&2
    """
}
