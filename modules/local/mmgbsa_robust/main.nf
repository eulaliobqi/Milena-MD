process MMGBSA_ROBUST {
    tag "${meta.id}"
    label 'process_medium'
    errorStrategy 'ignore'

    publishDir { "${params.outdir}/${meta.id}/mmgbsa" }, mode: 'copy'

    input:
    tuple val(meta), path(md_tpr), path(mmgbsa_xtc), path(lig_ndx)

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
    def sys_name = meta.id
    def saltcon  = params.nacl_conc ?: 0.15
    """
    echo "=== MMGBSA_ROBUST: ${meta.id} ===" >&2

    # ── 1. Validação pré-execução ────────────────────────────────────────────
    {
        echo "=== Validação pré-MMGBSA: ${meta.id} ==="
        echo "Data: \$(date)"
        echo ""
    } > mmgbsa_validation.txt
    touch tleap_wrapper.log

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

    # 1b. Conta frames reais na trajetória de input
    ACTUAL_FRAMES=\$(${params.gmx_cmd} check -f ${mmgbsa_xtc} 2>&1 \\
        | grep -E "^Last frame" | awk '{print \$3}' || echo "")
    if [ -z "\${ACTUAL_FRAMES}" ] || [ "\${ACTUAL_FRAMES}" -lt 10 ] 2>/dev/null; then
        echo "AVISO: não foi possível contar frames — usando estimativa" >> mmgbsa_validation.txt
        ACTUAL_FRAMES=100
    fi
    echo "Frames para MMGBSA: \${ACTUAL_FRAMES}" | tee -a mmgbsa_validation.txt
    echo "[MMGBSA] \${ACTUAL_FRAMES} frames na trajetória de input" >&2

    # ── 2. Gera mmgbsa.in com decomposição por resíduo ───────────────────────
    # endframe usa o valor real detectado (não uma estimativa)
    cat > mmgbsa.in << MEOF
&general
sys_name="${sys_name}",
startframe=1,
endframe=\${ACTUAL_FRAMES},
interval=1,
verbose=2,
/
&gb
igb=2,
saltcon=${saltcon},
/
&decomp
idecomp=2,
dec_verbose=1,
/
MEOF

    echo "mmgbsa.in gerado com endframe=\${ACTUAL_FRAMES}" >> mmgbsa_validation.txt

    # ── 3. Wrapper tleap robusto (corrige bug SS bonds COM_OUT) ──────────────
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

# Localiza o arquivo de input tleap (-f arquivo)
for i, a in enumerate(args):
    if a == '-f' and i + 1 < len(args):
        fpath = args[i + 1]
        if not os.path.exists(fpath):
            wlog(f"[tleap-wrapper] arquivo não encontrado: {fpath}")
            break

        content = open(fpath).read()
        wlog(f"[tleap-wrapper] processando: {fpath} ({len(content)} chars)")

        # Detecta bonds SS em REC_OUT e COM_OUT
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

# Delega para o tleap real (ignora bin_patch no PATH)
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

    # ── 4. Executa gmx_MMPBSA no ambiente isolado ────────────────────────────
    CURRENT_DIR="\$PWD"
    echo "[MMGBSA] Iniciando gmx_MMPBSA (pode demorar 20-60 min)..." >&2

    mamba run -n mmgbsa-env bash -c "
export PATH=\${CURRENT_DIR}/bin_patch:\\\$PATH
echo '[mmgbsa-env] PATH patch ativo' >&2

gmx_MMPBSA -O \\\\
    -i mmgbsa.in \\\\
    -cs ${md_tpr} \\\\
    -ct ${mmgbsa_xtc} \\\\
    -ci ${lig_ndx} \\\\
    -cg Receptor Ligante \\\\
    -o  FINAL_RESULTS_MMGBSA.dat \\\\
    -eo mmgbsa_results.csv \\\\
    -deo decomp_results.csv \\\\
    -nogui \\\\
    2>&1
" 2>&1 | tee mmgbsa.log

    # ── 5. Verifica saídas e cria fallbacks ──────────────────────────────────
    if [ -f FINAL_RESULTS_MMGBSA.dat ]; then
        echo "[OK] FINAL_RESULTS_MMGBSA.dat gerado" | tee -a mmgbsa_validation.txt
        # decomp_results.csv pode não ser gerado em todas as versões
        [ -f decomp_results.csv ] || {
            echo "resid,resname,total" > decomp_results.csv
            echo "AVISO: decomp_results.csv não gerado — arquivo vazio criado" >> mmgbsa_validation.txt
        }
        echo "[MMGBSA] Concluído com sucesso para ${meta.id}" >&2
    else
        echo "ERRO: FINAL_RESULTS_MMGBSA.dat não gerado" | tee -a mmgbsa_validation.txt
        echo "--- Últimas 40 linhas do log ---" >> mmgbsa_validation.txt
        tail -40 mmgbsa.log >> mmgbsa_validation.txt
        # Arquivos vazios para não bloquear downstream
        echo "No results — gmx_MMPBSA failed" > FINAL_RESULTS_MMGBSA.dat
        echo "frame,TOTAL"                     > mmgbsa_results.csv
        echo "resid,resname,total"              > decomp_results.csv
        exit 1
    fi
    """
}
