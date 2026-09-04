// MM-GBSA via gmx_MMPBSA. Este módulo nunca produziu um resultado real em
// nenhuma das 4 rodadas anteriores do projeto (DN773/DN1441/DN1937/DN2954 —
// todos os FINAL_RESULTS_MMGBSA.dat publicados são o placeholder de falha
// "No results", não números). Causa raiz e fix abaixo, confirmados por
// inspeção cruzada com o projeto irmão tatiane-MD (mesmo bug, mesmo fix,
// ver git log de modules/local/mmgbsa/main.nf lá: commits 2026-07-28):
//
// 1. INDIREÇÃO DUPLA DE CONDA (causa raiz de -cs/-ct/-ci chegarem vazios em
//    toda tentativa anterior): a task rodava no ambiente default e tentava
//    "alcançar" mmgbsa-env via `mamba run -n mmgbsa-env` DE DENTRO do
//    script — ativação aninhada sobre a ativação que o próprio Nextflow já
//    faz (process.conda), cenário onde $@ do subprocesso intermediário se
//    perde. Fix: `withName: 'MMGBSA_ROBUST' { conda = mmgbsa-env }` em
//    nextflow.config — Nextflow ativa o env nativamente, gmx_MMPBSA chamado
//    DIRETO abaixo, sem wrapper/mamba run.
// 2. Faltava -cp (topologia do complexo): gmx_MMPBSA em modo GROMACS exige
//    -cs (tpr) E -cp (o .top texto, não só o binário) — o .tpr sozinho não
//    basta para o parser ParmEd/GROMACS do gmx_MMPBSA.
// 3. -cg usava NOMES de grupo ("Receptor Ligante") — gmx_MMPBSA quer os
//    ÍNDICES NUMÉRICOS dos grupos no .ndx, na ordem em que aparecem
//    (convenção GROMACS 0..N), não os nomes.
process MMGBSA_ROBUST {
    tag "${meta.id}"
    label 'process_medium'
    errorStrategy 'ignore'

    publishDir { "${params.outdir}/${meta.id}/mmgbsa" }, mode: 'copy'

    input:
    tuple val(meta), path(md_tpr), path(mmgbsa_xtc), path(lig_ndx),
                     path(top, stageAs: 'input.top'), path(itps, stageAs: 'itp_in/*')

    output:
    tuple val(meta),
        path("mmgbsa_results.csv"),
        path("FINAL_RESULTS_MMGBSA.dat"),
        path("decomp_results.csv"),       emit: results
    tuple val(meta),
        path("mmgbsa.log"),
        path("mmgbsa_validation.txt"),
        path("tleap_wrapper.log"),        emit: logs

    script:
    def sys_name  = meta.id
    def saltcon   = params.nacl_conc ?: 0.15
    // CHARMM36/CGenFF não é topologia AMBER nativa -- gmx_MMPBSA só aceita
    // PBRadii=7 (charmm_radii) pra esse caso (GMXMMPBSA/make_top.py). Aqui
    // cobre 100% do sistema (receptor + peptídeo GORE1-2T são só os 20
    // aminoácidos padrão, charmm_radii define raio por nome de resíduo pra
    // todos eles) -- diferente do projeto irmão tatiane-MD, que precisou de
    // PBRadii=3 (mbondi2, por elemento) porque o ligante dela é CGenFF
    // (molécula pequena não-padrão, sem entrada em charmm_radii). Pra AMBER
    // (padrão deste pipeline) PBRadii fica de fora — gmx_MMPBSA usa o
    // default correto pra topologia amber nativa sem precisar da flag.
    def is_charmm    = params.forcefield.toString().startsWith('charmm36')
    def pbradii_line = is_charmm ? "PBRadii=7,\n" : ""
    """
    echo "=== MMGBSA_ROBUST: ${meta.id} ===" >&2

    {
        echo "=== Validação pré-MMGBSA: ${meta.id} ==="
        echo "Data: \$(date)"
        echo ""
    } > mmgbsa_validation.txt
    touch tleap_wrapper.log

    cp ${top} topol.top
    cp itp_in/*.itp .

    # 1a. Verifica grupos no lig.ndx (necessários para -cg Receptor Ligante)
    echo "--- Grupos em lig.ndx ---" >> mmgbsa_validation.txt
    echo q | ${params.gmx_cmd} make_ndx \\
        -f ${md_tpr} -n ${lig_ndx} \\
        2>&1 | grep -E "^ *[0-9]+" >> mmgbsa_validation.txt || true

    if ! grep -q "Receptor" ${lig_ndx}; then
        echo "ERRO FATAL: grupo 'Receptor' não encontrado em ${lig_ndx}" | tee -a mmgbsa_validation.txt
        exit 1
    fi
    if ! grep -q "Ligante" ${lig_ndx}; then
        echo "ERRO FATAL: grupo 'Ligante' não encontrado em ${lig_ndx}" | tee -a mmgbsa_validation.txt
        exit 1
    fi
    echo "OK: grupos 'Receptor' e 'Ligante' confirmados" >> mmgbsa_validation.txt

    # 1b. Índice NUMÉRICO dos grupos (gmx_MMPBSA quer "-cg <idx_rec> <idx_lig>",
    # não os nomes -- convenção GROMACS, blocos numerados 0..N na ordem do .ndx).
    IDX=0
    REC_IDX=""
    LIG_IDX=""
    while IFS= read -r line; do
        case "\$line" in
            \\[*\\])
                name=\$(echo "\$line" | sed -E 's/^\\[ *//; s/ *\\]\$//')
                [ "\$name" = "Receptor" ] && REC_IDX=\$IDX
                [ "\$name" = "Ligante" ]  && LIG_IDX=\$IDX
                IDX=\$((IDX + 1))
                ;;
        esac
    done < ${lig_ndx}
    echo "Índices calculados: Receptor=\${REC_IDX} Ligante=\${LIG_IDX}" | tee -a mmgbsa_validation.txt
    if [ -z "\${REC_IDX}" ] || [ -z "\${LIG_IDX}" ]; then
        echo "ERRO FATAL: não foi possível determinar índice numérico de Receptor/Ligante" | tee -a mmgbsa_validation.txt
        exit 1
    fi

    # 1c. Conta frames reais na trajetória de input
    ACTUAL_FRAMES=\$(${params.gmx_cmd} check -f ${mmgbsa_xtc} 2>&1 \\
        | grep -E "^Last frame" | awk '{print \$3}' || echo "")
    if [ -z "\${ACTUAL_FRAMES}" ] || [ "\${ACTUAL_FRAMES}" -lt 10 ] 2>/dev/null; then
        echo "AVISO: não foi possível contar frames — usando estimativa" >> mmgbsa_validation.txt
        ACTUAL_FRAMES=100
    fi
    echo "Frames para MMGBSA: \${ACTUAL_FRAMES}" | tee -a mmgbsa_validation.txt
    echo "[MMGBSA] \${ACTUAL_FRAMES} frames na trajetória de input" >&2

    # ── 2. Gera mmgbsa.in com decomposição por resíduo ───────────────────────
    # print_res="within 6": limita a decomposição aos resíduos na interface
    # (6 A do ligante) -- sem isso, idecomp=2 decompõe TODO o receptor
    # (235 resíduos), inviabilizando o tempo de execução sem necessidade.
    cat > mmgbsa.in << MEOF
&general
sys_name="${sys_name}",
startframe=1,
endframe=\${ACTUAL_FRAMES},
interval=1,
verbose=2,
${pbradii_line}/
&gb
igb=2,
saltcon=${saltcon},
/
&decomp
idecomp=2,
dec_verbose=1,
print_res="within 6",
/
MEOF
    echo "mmgbsa.in gerado com endframe=\${ACTUAL_FRAMES}" >> mmgbsa_validation.txt

    # ── 3. Wrapper tleap robusto (corrige bug SS bonds COM_OUT) ──────────────
    # Real pra este receptor (arquitetura de dissulfetos conhecida da
    # tripsina, 6 CYS/sistema) -- mantido de versões anteriores deste módulo,
    # nunca chegou a ser exercido em produção porque o bug #1 (indireção de
    # conda) sempre abortava antes de tleap rodar.
    mkdir -p bin_patch
    cat > bin_patch/tleap << 'WEOF'
#!/usr/bin/env python3
# Wrapper tleap para gmx_MMPBSA 1.6.x
# Corrige bug: indices SS bonds em COM_OUT gerados com offset errado.
# O fix copia os indices corretos de REC_OUT para COM_OUT.
# Registra tudo em tleap_wrapper.log para diagnostico.
import sys, os, re, subprocess

LOG = os.path.join(os.getcwd(), 'tleap_wrapper.log')

def wlog(msg):
    with open(LOG, 'a') as f:
        f.write(msg + '\n')

args = sys.argv[1:]
wlog(f"[tleap-wrapper] chamado com: {' '.join(args)}")

for i, a in enumerate(args):
    if a == '-f' and i + 1 < len(args):
        fpath = args[i + 1]
        if not os.path.exists(fpath):
            wlog(f"[tleap-wrapper] arquivo não encontrado: {fpath}")
            break

        content = open(fpath).read()
        wlog(f"[tleap-wrapper] processando: {fpath} ({len(content)} chars)")

        rec_bonds = re.findall(
            r'bond REC_OUT\.(\d+)\.SG\s+REC_OUT\.(\d+)\.SG', content)
        com_bonds = re.findall(
            r'bond COM_OUT\.(\d+)\.SG\s+COM_OUT\.(\d+)\.SG', content)

        wlog(f"[tleap-wrapper] REC SS bonds: {rec_bonds}")
        wlog(f"[tleap-wrapper] COM SS bonds: {com_bonds}")

        if rec_bonds and com_bonds:
            if len(com_bonds) != len(rec_bonds):
                wlog(f"[tleap-wrapper] AVISO: {len(com_bonds)} COM vs "
                     f"{len(rec_bonds)} REC bonds — aplicando zip truncado")

            modified = content
            fixes = 0
            for (cw0, cw1), (rr0, rr1) in zip(com_bonds, rec_bonds):
                old = f'bond COM_OUT.{cw0}.SG COM_OUT.{cw1}.SG'
                new = f'bond COM_OUT.{rr0}.SG COM_OUT.{rr1}.SG'
                if old != new:
                    modified = modified.replace(old, new, 1)
                    fixes += 1
                    wlog(f"[tleap-wrapper] FIX {fixes}: '{old}' → '{new}'")

            if fixes > 0:
                open(fpath, 'w').write(modified)
                wlog(f"[tleap-wrapper] {fixes} correção(ões) aplicada(s)")
            else:
                wlog("[tleap-wrapper] índices já corretos, sem modificação")
        elif not rec_bonds and not com_bonds:
            wlog("[tleap-wrapper] nenhuma ponte SS detectada — sem correção necessária")

        break

path_dirs = os.environ.get('PATH', '').split(':')
for d in path_dirs:
    if 'bin_patch' in d:
        continue
    for name in ('tleap', 'teLeap'):
        exe = os.path.join(d, name)
        if os.path.isfile(exe) and os.access(exe, os.X_OK):
            wlog(f"[tleap-wrapper] executando tleap real: {exe}")
            ret = subprocess.run([exe] + args)
            wlog(f"[tleap-wrapper] retcode: {ret.returncode}")
            sys.exit(ret.returncode)

wlog("[tleap-wrapper] ERRO FATAL: nenhum tleap real encontrado no PATH")
sys.exit(1)
WEOF
    chmod +x bin_patch/tleap
    echo "bin_patch/tleap criado" >> mmgbsa_validation.txt

    # ── 4. Executa gmx_MMPBSA — chamada DIRETA (Nextflow já ativou
    # mmgbsa-env nativamente via process.conda para esta task, ver
    # nextflow.config), sem mamba run/wrapper externo. -cp topol.top e -cg
    # com índices numéricos são os dois fixes reais desta rodada (ver
    # comentário do módulo).
    export PATH="\$PWD/bin_patch:\$PATH"
    echo "[MMGBSA] Iniciando gmx_MMPBSA (pode demorar 20-60 min)..." >&2

    gmx_MMPBSA -O \\
        -i   mmgbsa.in \\
        -cs  ${md_tpr} \\
        -ct  ${mmgbsa_xtc} \\
        -ci  ${lig_ndx} \\
        -cg  \${REC_IDX} \${LIG_IDX} \\
        -cp  topol.top \\
        -o   FINAL_RESULTS_MMGBSA.dat \\
        -eo  mmgbsa_results.csv \\
        -deo decomp_results.csv \\
        -nogui \\
        2>&1 | tee mmgbsa.log

    # ── 5. Verifica saídas e cria fallbacks ──────────────────────────────────
    if [ -s FINAL_RESULTS_MMGBSA.dat ] && ! grep -q "^No results" FINAL_RESULTS_MMGBSA.dat; then
        echo "[OK] FINAL_RESULTS_MMGBSA.dat gerado" | tee -a mmgbsa_validation.txt
        [ -f decomp_results.csv ] || {
            echo "resid,resname,total" > decomp_results.csv
            echo "AVISO: decomp_results.csv não gerado — arquivo vazio criado" >> mmgbsa_validation.txt
        }
        echo "[MMGBSA] Concluído com sucesso para ${meta.id}" >&2
    else
        echo "ERRO: FINAL_RESULTS_MMGBSA.dat não gerado" | tee -a mmgbsa_validation.txt
        echo "--- Últimas 40 linhas do log ---" >> mmgbsa_validation.txt
        tail -40 mmgbsa.log >> mmgbsa_validation.txt
        # Arquivos vazios para não bloquear downstream — sair com 0 é
        # proposital: com errorStrategy 'ignore', um exit != 0 faz a tarefa
        # não emitir NENHUM output (mesmo os placeholders acima existindo no
        # work dir), o que travava PLOT rio abaixo. O log completo continua
        # em mmgbsa_validation.txt/mmgbsa.log para diagnóstico.
        echo "No results — gmx_MMPBSA failed" > FINAL_RESULTS_MMGBSA.dat
        echo "frame,TOTAL"                     > mmgbsa_results.csv
        echo "resid,resname,total"              > decomp_results.csv
    fi
    """
}
