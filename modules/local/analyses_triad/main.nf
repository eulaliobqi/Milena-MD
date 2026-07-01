process ANALYSES_TRIAD {
    tag "${meta.id}"
    label 'process_low'

    publishDir { "${params.outdir}/${meta.id}/analise" }, mode: 'copy'

    input:
    tuple val(meta), path(md_tpr), path(md_fit_xtc), path(lig_ndx)

    output:
    tuple val(meta), path("dist_r1.xvg"), path("dist_r2.xvg"),
                     path("dist_r3.xvg"), path("dist_r4.xvg"),
                     path("sasa_r1.xvg"), path("sasa_r2.xvg"),
                     path("sasa_r3.xvg"), path("sasa_r4.xvg"), emit: triad
    tuple val(meta), path("triad_info.txt"), emit: info

    script:
    def r1 = meta.triad_1
    def r2 = meta.triad_2
    def r3 = meta.triad_3
    def r4 = meta.triad_4
    """
    echo "=== ANALYSES_TRIAD: ${meta.id} ===" >&2
    echo "Resíduos de interesse: ${r1} ${r2} ${r3} ${r4}" >&2

    # Grava rótulos para plot_results.py
    printf '${r1}\\n${r2}\\n${r3}\\n${r4}\\n' > triad_info.txt

    # Conta grupos existentes no lig.ndx para determinar próximos índices
    N_CURR=\$(echo q | ${params.gmx_cmd} make_ndx \\
        -f ${md_tpr} -n ${lig_ndx} -o _tmp_count.ndx 2>&1 \\
        | grep -cE "^ *[0-9]+ +[A-Za-z]")
    rm -f _tmp_count.ndx
    R1_IDX=\${N_CURR}
    R2_IDX=\$((N_CURR + 1))
    R3_IDX=\$((N_CURR + 2))
    R4_IDX=\$((N_CURR + 3))

    # Adiciona os 4 resíduos ao ndx
    ${params.gmx_cmd} make_ndx -f ${md_tpr} -n ${lig_ndx} -o triad.ndx << MNDX
r ${r1}
name \${R1_IDX} Res1_Cat
r ${r2}
name \${R2_IDX} Res2_Cat
r ${r3}
name \${R3_IDX} Res3_Cat
r ${r4}
name \${R4_IDX} Res4_Cat
q
MNDX

    printf 'Ligante\\nRes1_Cat\\n' | ${params.gmx_cmd} mindist \\
        -s ${md_tpr} -f ${md_fit_xtc} -n triad.ndx -od dist_r1.xvg -tu ns

    printf 'Ligante\\nRes2_Cat\\n' | ${params.gmx_cmd} mindist \\
        -s ${md_tpr} -f ${md_fit_xtc} -n triad.ndx -od dist_r2.xvg -tu ns

    printf 'Ligante\\nRes3_Cat\\n' | ${params.gmx_cmd} mindist \\
        -s ${md_tpr} -f ${md_fit_xtc} -n triad.ndx -od dist_r3.xvg -tu ns

    printf 'Ligante\\nRes4_Cat\\n' | ${params.gmx_cmd} mindist \\
        -s ${md_tpr} -f ${md_fit_xtc} -n triad.ndx -od dist_r4.xvg -tu ns

    # SASA por resíduo catalítico (surface = proteína completa; output = resíduo individual)
    # Valores baixos indicam resíduo enterrado/em contato com ligante
    printf 'Protein\\nRes1_Cat\\n' | ${params.gmx_cmd} sasa \\
        -s ${md_tpr} -f ${md_fit_xtc} -n triad.ndx -o sasa_r1.xvg -tu ns

    printf 'Protein\\nRes2_Cat\\n' | ${params.gmx_cmd} sasa \\
        -s ${md_tpr} -f ${md_fit_xtc} -n triad.ndx -o sasa_r2.xvg -tu ns

    printf 'Protein\\nRes3_Cat\\n' | ${params.gmx_cmd} sasa \\
        -s ${md_tpr} -f ${md_fit_xtc} -n triad.ndx -o sasa_r3.xvg -tu ns

    printf 'Protein\\nRes4_Cat\\n' | ${params.gmx_cmd} sasa \\
        -s ${md_tpr} -f ${md_fit_xtc} -n triad.ndx -o sasa_r4.xvg -tu ns

    echo "[OK] Distâncias e SASA catalíticos concluídos para ${meta.id}" >&2
    """
}
