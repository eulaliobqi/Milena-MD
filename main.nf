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
    // tríade + MM-GBSA). O CSV bruto do MM-GBSA (mmgbsa_results.csv) vem de
    // MMGBSA_ROBUST.out.results — é o único dos dois módulos que emite o CSV
    // que plot_results.py consome via --mmgbsa-csv (MMGBSA_INTERPRET só emite
    // painel_mmgbsa.png + mmgbsa_summary.txt, não um CSV).
    //
    // STABILITY_FILTER/CLUSTERING/MMGBSA_ROBUST têm errorStrategy 'ignore'
    // (gmx_MMPBSA é sabidamente instável) — se qualquer um falhar, essa amostra
    // não emite em MMGBSA_ROBUST.out.results. Para não travar o PLOT principal
    // (RMSD/RMSF/Rg/tríade, que já são análises validadas) nesse cenário,
    // usamos join(remainder: true) ancorado no canal de amostras de entrada e
    // caímos para um CSV placeholder vazio (frame,TOTAL) quando o MM-GBSA não
    // produziu resultado — plot_results.py trata isso como "sem MM-GBSA" e
    // simplesmente omite o painel de energia de ligação.
    ch_sample_ids = ch_input.map { meta, receptor, ligand -> meta }

    ch_mmgbsa_csv_for_plot = ch_sample_ids
        .join(MMGBSA_ROBUST.out.results.map { meta, csv, dat, decomp -> tuple(meta, csv) }, remainder: true)
        .map { meta, csv -> tuple(meta, csv ?: file("${projectDir}/assets/no_mmgbsa.csv")) }

    ch_plot_xvg = ANALYSES.out.xvg
        .join(ANALYSES_TRIAD.out.triad, by: [0])
        .join(ANALYSES_TRIAD.out.info,  by: [0])
        .map { meta, xvgs, ndx, d1, d2, d3, d4, s1, s2, s3, s4, info ->
            def all_files = (xvgs instanceof List ? xvgs : [xvgs]) +
                            [d1, d2, d3, d4, s1, s2, s3, s4, info]
            tuple(meta, all_files, ndx)
        }

    // NOTA: meta em ch_plot_xvg vem de ANALYSES_TRIAD (contém triad_1..4),
    // enquanto meta em ch_mmgbsa_csv_for_plot vem do canal base (só "id") —
    // são mapas diferentes em conteúdo. Para não depender de igualdade de Map
    // (frágil e não garantida pelo operador join por valor de mapa completo),
    // o merge final usa meta.id (String) como chave explícita.
    ch_plot_xvg_keyed   = ch_plot_xvg.map           { meta, all_files, ndx -> tuple(meta.id, meta, all_files, ndx) }
    ch_mmgbsa_csv_keyed = ch_mmgbsa_csv_for_plot.map { meta, csv           -> tuple(meta.id, csv) }

    ch_plot_input = ch_plot_xvg_keyed
        .join(ch_mmgbsa_csv_keyed, by: 0)
        .map { id, meta, all_files, ndx, csv -> tuple(meta, all_files, ndx, csv) }

    PLOT(ch_plot_input)
}
