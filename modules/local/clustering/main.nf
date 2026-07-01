process CLUSTERING {
    tag "${meta.id}"
    label 'process_medium'
    errorStrategy 'ignore'

    publishDir { "${params.outdir}/${meta.id}/mmgbsa_prep" }, mode: 'copy'

    input:
    tuple val(meta), path(md_tpr), path(stable_xtc), path(lig_ndx)

    output:
    tuple val(meta), path(md_tpr), path("mmgbsa_input.xtc"), path(lig_ndx), emit: for_mmgbsa
    tuple val(meta), path("clusterid.xvg"), path("cluster_centers.pdb"),    emit: clusters
    tuple val(meta), path("clustering_report.txt"),                          emit: report

    script:
    def cutoff        = params.cluster_cutoff             ?: 0.2
    def n_cl          = params.mmgbsa_n_clusters          ?: 3
    def fpc           = params.mmgbsa_frames_per_cluster  ?: 50
    def target_frames = n_cl * fpc
    """
    echo "=== CLUSTERING: ${meta.id} ===" >&2

    {
        echo "=== Relatório de Clustering ==="
        echo "Sistema  : ${meta.id}"
        echo "Método   : GROMOS"
        echo "Grupo    : Ligante (pose de ligação)"
        echo "Cutoff   : ${cutoff} nm"
        echo ""
    } > clustering_report.txt

    # ── Clustering GROMOS pelo RMSD do ligante ─────────────────────────────
    # Usar o ligante como referência captura variações da pose de ligação —
    # mais informativo biologicamente do que o backbone completo.
    printf 'Ligante\\nLigante\\n' | ${params.gmx_cmd} cluster \\
        -s ${md_tpr} \\
        -f ${stable_xtc} \\
        -n ${lig_ndx} \\
        -method gromos \\
        -cutoff ${cutoff} \\
        -o clusters.xpm \\
        -g cluster.log \\
        -clid clusterid.xvg \\
        -cl cluster_centers.pdb \\
        -tu ns \\
        2>&1 | tee -a clustering_report.txt

    # Extrai estatísticas do cluster.log
    echo "" >> clustering_report.txt
    echo "--- Clusters detectados ---" >> clustering_report.txt
    grep -E "^(cl\\.|Total|Middle|Found)" cluster.log 2>/dev/null \\
        | head -50 >> clustering_report.txt || true

    # ── Cria subtrajetória para MMGBSA por subsampling da fase estável ───────
    # Conta frames reais na fase estável
    STABLE_FRAMES=\$(${params.gmx_cmd} check -f ${stable_xtc} 2>&1 \\
        | grep -E "^Last frame" | awk '{print \$3}' || echo "1000")

    # Calcula stride para atingir target_frames
    SKIP=\$(python3 -c "print(max(1, int(\${STABLE_FRAMES} // ${target_frames})))")

    echo "" >> clustering_report.txt
    echo "--- Subsampling para MMGBSA ---" >> clustering_report.txt
    echo "Frames estáveis : \${STABLE_FRAMES}" >> clustering_report.txt
    echo "Target frames   : ${target_frames}  (${n_cl} clusters × ${fpc} frames)" >> clustering_report.txt
    echo "Stride usado    : \${SKIP}" >> clustering_report.txt

    EST_FRAMES=\$(( \${STABLE_FRAMES} / \${SKIP} ))
    echo "[CLUSTER] stride=\${SKIP} → ~\${EST_FRAMES} frames para MMGBSA" >&2

    echo "System" | ${params.gmx_cmd} trjconv \\
        -s ${md_tpr} \\
        -f ${stable_xtc} \\
        -o mmgbsa_input.xtc \\
        -skip \${SKIP} \\
        2>&1 | tail -5

    # Verifica output
    N_OUT=\$(${params.gmx_cmd} check -f mmgbsa_input.xtc 2>&1 \\
        | grep -E "^Last frame" | awk '{print \$3}' || echo "?")
    echo "Frames para MMGBSA : \${N_OUT}" | tee -a clustering_report.txt
    echo "[OK] mmgbsa_input.xtc: \${N_OUT} frames" >&2
    """
}
