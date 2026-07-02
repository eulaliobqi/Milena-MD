#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { PREPARE_PH        } from './modules/local/prepare_ph/main.nf'
include { PREPARE_COMPLEX   } from './modules/local/prepare_complex/main.nf'
include { TOPOLOGY          } from './modules/local/topology/main.nf'
include { BOX_SOLVATE_IONS  } from './modules/local/box_solvate_ions/main.nf'
include { MINIMIZATION      } from './modules/local/minimization/main.nf'
include { NVT               } from './modules/local/nvt/main.nf'
include { NPT               } from './modules/local/npt/main.nf'
include { PRODUCTION        } from './modules/local/production/main.nf'
include { POSTPROCESS       } from './modules/local/postprocess/main.nf'
include { ANALYSES          } from './modules/local/analyses/main.nf'
include { ANALYSES_TRIAD    } from './modules/local/analyses_triad/main.nf'
include { STABILITY_FILTER  } from './modules/local/stability_filter/main.nf'
include { CLUSTERING        } from './modules/local/clustering/main.nf'
include { MMGBSA_ROBUST     } from './modules/local/mmgbsa_robust/main.nf'
include { MMGBSA_INTERPRET  } from './modules/local/mmgbsa_interpret/main.nf'
include { PLOT              } from './modules/local/plot/main.nf'

workflow {
    if (!params.input) {
        error "Informe o samplesheet: --input samplesheet.csv"
    }

    // Canal principal — meta sem campos de tríade para preservar cache upstream
    Channel
        .fromPath(params.input, checkIfExists: true)
        .splitCsv(header: true, strip: true)
        .map { row ->
            def meta     = [id: row.sample_id]
            def receptor = file(row.receptor, checkIfExists: true)
            def ligand   = file(row.ligand,   checkIfExists: true)
            tuple(meta, receptor, ligand)
        }
        .set { ch_input }

    // Canal de resíduos de interesse — injetado apenas em ANALYSES_TRIAD
    Channel
        .fromPath(params.input, checkIfExists: true)
        .splitCsv(header: true, strip: true)
        .map { row ->
            tuple(
                [id: row.sample_id],
                row.triad_1 ?: params.triad_1,
                row.triad_2 ?: params.triad_2,
                row.triad_3 ?: params.triad_3,
                row.triad_4 ?: params.triad_4
            )
        }
        .set { ch_triad_params }

    // ── Cadeia principal de preparação + MD ────────────────────────────────────
    PREPARE_PH(ch_input)
    PREPARE_COMPLEX(PREPARE_PH.out.protonated)
    TOPOLOGY(PREPARE_COMPLEX.out.complexo)
    BOX_SOLVATE_IONS(TOPOLOGY.out.topology)
    MINIMIZATION(BOX_SOLVATE_IONS.out.system)
    NVT(MINIMIZATION.out.system)
    NPT(NVT.out.system)
    PRODUCTION(NPT.out.system)
    POSTPROCESS(PRODUCTION.out.traj)

    // Junta complexo.pdb (preparação) com trajetória pós-processada
    ch_analyses = PREPARE_COMPLEX.out.complexo
        .join(POSTPROCESS.out.fit, by: [0])

    ANALYSES(ch_analyses)

    // ── Análises estendidas (distâncias tríade + bolsão S1) ───────────────────
    ch_extended = POSTPROCESS.out.fit
        .join(ANALYSES.out.xvg, by: [0])
        .map { meta, tpr, xtc, xvgs, ndx -> tuple(meta, tpr, xtc, ndx) }

    // Injeta resíduos de interesse no canal só para ANALYSES_TRIAD
    ch_triad_input = ch_extended
        .join(ch_triad_params, by: [0])
        .map { meta, tpr, xtc, ndx, t1, t2, t3, t4 ->
            tuple(meta + [triad_1: t1, triad_2: t2, triad_3: t3, triad_4: t4], tpr, xtc, ndx)
        }

    ANALYSES_TRIAD(ch_triad_input)

    // ══════════════════════════════════════════════════════════════════════════
    // Pipeline MM-GBSA robusto (STABILITY_FILTER → CLUSTERING → MMGBSA_ROBUST →
    // MMGBSA_INTERPRET), encadeado automaticamente no fluxo principal — ao
    // contrário do MD-gromacs original, onde estes 4 módulos existem mas só
    // rodam manualmente via scripts/run_mmgbsa.sh.
    // ══════════════════════════════════════════════════════════════════════════

    // STABILITY_FILTER espera: meta, md_tpr, md_xtc, xvg_files, lig_ndx
    // - md_tpr/md_xtc: trajetória alinhada pós-processada (POSTPROCESS.out.fit
    //   emite md.tpr + md_fit.xtc — reaproveitados aqui como tpr/xtc de entrada)
    // - xvg_files/lig_ndx: saída de ANALYSES (contém rmsd_backbone.xvg, usado
    //   pelo stability_filter.py para detectar o platô de RMSD)
    ch_stability_input = POSTPROCESS.out.fit
        .join(ANALYSES.out.xvg, by: [0])
        .map { meta, tpr, xtc, xvgs, ndx -> tuple(meta, tpr, xtc, xvgs, ndx) }

    STABILITY_FILTER(ch_stability_input)

    // CLUSTERING espera: meta, md_tpr, stable_xtc, lig_ndx
    CLUSTERING(STABILITY_FILTER.out.stable)

    // MMGBSA_ROBUST espera: meta, md_tpr, mmgbsa_xtc, lig_ndx
    MMGBSA_ROBUST(CLUSTERING.out.for_mmgbsa)

    // MMGBSA_INTERPRET espera: meta, mmgbsa_csv, mmgbsa_dat, decomp_csv, stability_report
    ch_mmgbsa_interpret_input = MMGBSA_ROBUST.out.results
        .join(STABILITY_FILTER.out.report, by: [0])
        .map { meta, csv, dat, decomp, report -> tuple(meta, csv, dat, decomp, report) }

    MMGBSA_INTERPRET(ch_mmgbsa_interpret_input)

    // ── Plot final — painel composto (RMSD/RMSF/Rg/contatos/H-bonds/SASA +
    // tríade [+ MM-GBSA se disponível]).
    //
    // PLOT NÃO depende de MMGBSA_ROBUST via canal/join — gmx_MMPBSA se
    // mostrou instável nas primeiras rodadas de DN773-GORE12T (2026-07-02:
    // já falhou por args vazios, args duplicados e erro de ambiente), e um
    // join anterior com MMGBSA_ROBUST.out.results (mesmo com
    // remainder:true + CSV placeholder) travou o PLOT inteiro quando
    // MM-GBSA falhava — ver commit 52f605d. Em vez de religar esse join,
    // PLOT usa sempre assets/no_mmgbsa.csv (arquivo estático, sem canal);
    // plot_results.py já trata esse CSV placeholder graciosamente
    // (load_mmgbsa_csv retorna None com <2 linhas → has_mmgbsa=False, painel
    // MM-GBSA simplesmente omitido). Resultado: RMSD/RMSF/Rg/tríade sempre
    // saem, MM-GBSA é estritamente bônus quando MMGBSA_INTERPRET funcionar
    // (painel separado, publicado em mmgbsa/painel_mmgbsa.png).
    ch_no_mmgbsa = file("${projectDir}/assets/no_mmgbsa.csv")

    ch_plot_input = ANALYSES.out.xvg
        .join(ANALYSES_TRIAD.out.triad, by: [0])
        .join(ANALYSES_TRIAD.out.info,  by: [0])
        .map { meta, xvgs, ndx, d1, d2, d3, d4, s1, s2, s3, s4, info ->
            def all_files = (xvgs instanceof List ? xvgs : [xvgs]) +
                            [d1, d2, d3, d4, s1, s2, s3, s4, info]
            tuple(meta, all_files, ndx, ch_no_mmgbsa)
        }

    PLOT(ch_plot_input)
}
