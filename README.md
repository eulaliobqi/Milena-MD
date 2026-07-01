# Milena-MD — DN773-GORE12T (Nextflow DSL2 + MM-GBSA automatizado)

Pipeline Nextflow DSL2 para dinâmica molecular do complexo **DN773 (receptor,
tripsina digestiva de *Spodoptera*) + GORE12T (peptídeo/proteína ligante, 75
resíduos)**. É a 6ª isoforma de receptor da série iniciada em
[MD-gromacs](https://github.com/eulaliobqi/MD-gromacs) (ACR157 ✅, QCL936 ✅,
XP273 ✅, XP352 ❌ dissociou), e é um projeto **autônomo**: copia/adapta os
módulos já validados do MD-gromacs (não depende dele em runtime).

Principal diferença em relação ao MD-gromacs: aqui o **MM-GBSA robusto**
(`gmx_MMPBSA`, modelo GB igb=2, decomposição por resíduo idecomp=2) está
**encadeado no fluxo automático principal** (`main.nf`), em vez de rodar
manualmente via `scripts/run_mmgbsa.sh` como no MD-gromacs original.

---

## 1. Visão geral

- **Receptor:** DN773, cadeia A, resíduos 24–257
- **Ligante:** GORE12T, cadeia B, resíduos 1–75 (peptídeo/proteína — via
  caminho proteína-proteína normal, **sem** ACPYPE/small-molecule)
- **Tríade catalítica (DN773):** His69, Asp116, Ser213, Gly207 (bolsão S1)
- **Condições físicas:** pH 8,2 (pdb2pqr + PROPKA), 300 K, KCl 0,10 M,
  AMBER99SB-ILDN + TIP3P, caixa cúbica 2,0 nm
- **Produção:** 100 ns, Parrinello-Rahman, GPU (`gmx_mpi` + `mpirun -np 1`)
- **MM-GBSA:** filtragem de estabilidade (platô de RMSD) → clustering GROMOS
  (cutoff 0,2 nm, 3 clusters × 50 frames) → `gmx_MMPBSA` GB (igb=2,
  saltcon=0,10 M) + decomposição por resíduo → painel interpretativo

---

## 2. Estrutura do repositório

```
Milena-MD/
├── main.nf                              # Workflow — DAG completo (15 processos)
├── nextflow.config                      # Params + profiles (local/slurm/conda)
├── conf/{base,local,slurm}.config       # Recursos por label / executor
├── environment.yml                      # env md-gromacs (GROMACS CUDA, pdb2pqr...)
├── environment-mmgbsa.yml               # env mmgbsa-env isolado (AmberTools+gmx_MMPBSA)
│
├── modules/local/
│   ├── prepare_ph/          # pdb2pqr (PROPKA) — protonação receptor+ligante @ pH
│   ├── prepare_complex/     # Monta complexo.pdb (CYX, chains A/B)
│   ├── topology/            # gmx pdb2gmx → complexo.gro + topol.top + *.itp
│   ├── box_solvate_ions/    # editconf (cúbico 2 nm) + solvate + genion (KCl)
│   ├── minimization/        # EM steep
│   ├── nvt/                 # Equil. NVT 200 ps (V-rescale)
│   ├── npt/                 # Equil. NPT 500 ps (Berendsen)
│   ├── production/          # DM de produção 100 ns (Parrinello-Rahman)
│   ├── postprocess/         # Correção PBC + fit rot+trans
│   ├── analyses/            # RMSD bb/ligante, RMSF, Rg, contatos, H-bonds, SASA×2
│   ├── analyses_triad/      # Distâncias + SASA por resíduo catalítico (4)
│   ├── stability_filter/    # Detecta platô de RMSD → stable.xtc
│   ├── clustering/          # GROMOS + subamostragem → mmgbsa_input.xtc
│   ├── mmgbsa_robust/       # gmx_MMPBSA GB + decomposição (com wrapper tleap)
│   ├── mmgbsa_interpret/    # Painel + resumo interpretativo MM-GBSA
│   └── plot/                # Painel composto final (inclui --mmgbsa-csv)
│
├── bin/                                 # Scripts chamados pelos módulos acima
│   ├── pdb2pqr_process.py               # Pós-processa saída do pdb2pqr
│   ├── prepare_complex.py               # Monta complexo.pdb (CYX/chains)
│   ├── stability_filter.py              # Detecção de platô de RMSD
│   ├── mmgbsa_interpret.py              # Painel/resumo MM-GBSA
│   └── plot_results.py                  # Painel composto de figuras
│
├── data/
│   ├── DN773-receptor.pdb               # Cadeia A, resíduos 24-257
│   └── GORE12T-ligand.pdb               # Cadeia B, resíduos 1-75
│
├── assets/
│   └── samplesheet_dn773_gore12t.csv    # Samplesheet de entrada
│
└── complexo-DN773-GORE12T.pdb           # Pose de referência/preview (NÃO é input do pipeline)
```

---

## 3. Fluxo do pipeline (DAG completo)

```
PREPARE_PH → PREPARE_COMPLEX → TOPOLOGY → BOX_SOLVATE_IONS
→ MINIMIZATION → NVT (200 ps) → NPT (500 ps, Berendsen)
→ PRODUCTION (100 ns, Parrinello-Rahman) → POSTPROCESS
→ ANALYSES (RMSD bb+ligante, RMSF, Rg, contatos/mindist, H-bonds, SASA×2)
→ ANALYSES_TRIAD (dist + SASA por resíduo His69/Asp116/Ser213/Gly207)
→ STABILITY_FILTER (detecta platô de estabilidade do RMSD)
→ CLUSTERING (GROMOS, cutoff 0,2 nm, subamostragem 3×50 frames)
→ MMGBSA_ROBUST (gmx_MMPBSA GB, igb=2, saltcon=0,10, + decomposição idecomp=2)
→ MMGBSA_INTERPRET (painel + resumo)
→ PLOT (painel composto final, incluindo --mmgbsa-csv)
```

`complexo-DN773-GORE12T.pdb` na raiz é apenas uma pose de referência montada
manualmente para conferência visual — o pipeline recebe receptor e ligante
**separados** via samplesheet (é o que `PREPARE_PH`/`PREPARE_COMPLEX` esperam).

### Robustez do MM-GBSA

`STABILITY_FILTER`, `CLUSTERING` e `MMGBSA_ROBUST` têm `errorStrategy 'ignore'`
(gmx_MMPBSA é sabidamente instável — daí o wrapper `tleap` embutido em
`mmgbsa_robust` que corrige o bug de offset de índices de ponte dissulfeto na
gmx_MMPBSA 1.6.x). Se qualquer etapa do MM-GBSA falhar para uma amostra, o
`PLOT` dessa amostra simplesmente não roda (join direto entre `ANALYSES`/
`ANALYSES_TRIAD` e `MMGBSA_ROBUST.out.results`, sem fallback) — mesmo
comportamento padrão do resto da cadeia quando um processo falha.

---

## 4. Samplesheet

`assets/samplesheet_dn773_gore12t.csv`:

```csv
sample_id,receptor,ligand,triad_1,triad_2,triad_3,triad_4
dn773-gore12t,/home/eulalio/gromacs/Milena-MD/data/DN773-receptor.pdb,/home/eulalio/gromacs/Milena-MD/data/GORE12T-ligand.pdb,69,116,213,207
```

Colunas `triad_1..4` seguem a convenção His, Asp, Ser, S1 (mesma ordem usada
em todos os samplesheets do MD-gromacs). Caminhos devem ser **absolutos no
servidor**.

---

## 5. Execução no servidor

```bash
cd ~/gromacs/Milena-MD && git pull origin main
screen -S dn773-gore12t
nextflow run main.nf \
    --input  /home/eulalio/gromacs/Milena-MD/assets/samplesheet_dn773_gore12t.csv \
    --outdir /home/eulalio/gromacs/Milena-MD/DN773-GORE12T/MD \
    --pH 8.2 --time_ns 100 -profile local,conda
# Ctrl+A D para desanexar
```

Para reanexar depois: `screen -r dn773-gore12t`.

Ambientes necessários no servidor (criar antes, se ainda não existirem):

```bash
mamba env create -f environment.yml          # env md-gromacs
mamba env create -f environment-mmgbsa.yml   # env mmgbsa-env (isolado)
```

### Regras críticas (não-negociáveis)

| Regra | Por quê |
|---|---|
| `screen -S <nome>` ANTES do `nextflow run` | SIGTTOU mata o processo em background assim que a sessão SSH desconecta |
| `mamba`, nunca `conda` | Resolver de dependências mais rápido, ambientes mais estáveis |
| `gmx_mpi` + `mpirun -np 1 gmx_mpi ...` | Build CUDA-MPI do GROMACS 2026.0; `gmx` puro não existe no servidor |
| `python=3.11` no env `mmgbsa-env` | `gmx_MMPBSA`/AmberTools não têm wheels para Python 3.14 |
| `git pull` antes de rodar | Sempre sincronizar com o último fix commitado |
| `--titration-state-method propka` no pdb2pqr | pdb2pqr 3.x renomeou o argumento (era `--with-ph` sozinho em versões antigas) |
| `${System.getenv('HOME')}` no `nextflow.config` | Nextflow 26.x não resolve `$HOME` diretamente em `params`/`profiles` |
| `df -h /home` antes de nova simulação MD | Partição já chegou a 96% por causa de `af3_databases` — conferir espaço livre antes de rodar |

---

## 6. Verificação pós-execução

- **Preparação:** `<outdir>/dn773-gore12t/prep_ph/receptor_ph.pdb` e
  `ligand_ph.pdb` existem e têm átomos de H adicionados.
- **Topologia:** `<outdir>/dn773-gore12t/topo/topol.top` sem erros de
  `pdb2gmx.log` (conferir ausência de "Fatal error").
- **Produção:** `<outdir>/dn773-gore12t/prod/md.gro` existe e
  `gmx_mpi check -f md.xtc` reporta ~50.000 frames (100 ns / 2 ps de
  `nstxout-compressed`).
- **Análises:** `<outdir>/dn773-gore12t/analise/rmsd_backbone.xvg`,
  `rmsd_ligante.xvg`, `dist_r1..4.xvg` presentes e não vazios.
- **MM-GBSA:** `<outdir>/dn773-gore12t/mmgbsa/FINAL_RESULTS_MMGBSA.dat` com
  `ΔG bind (TOTAL)` reportado; `mmgbsa_summary.txt` com interpretação textual.
- **Painel final:** `<outdir>/dn773-gore12t/analise/painel_completo.png`
  (ou `painel_mmgbsa.png` em `mmgbsa/`) — abrir e conferir se o painel de
  energia de ligação foi incluído (se ausente, checar `mmgbsa.log` e
  `mmgbsa_validation.txt` para causa da falha).

Se precisar regenerar apenas os gráficos manualmente (útil se `PLOT` não
capturou o CSV de MM-GBSA por algum motivo):

```bash
mamba activate md-gromacs
python3 bin/plot_results.py \
    --analise-dir  <outdir>/dn773-gore12t/analise \
    --titulo       "DN773-GORE12T — DM 100 ns @ pH 8.2" \
    --window-ns    5 \
    --mmgbsa-csv   <outdir>/dn773-gore12t/mmgbsa/mmgbsa_results.csv
```

---

## 7. Próximos passos após a rodada

1. Conferir `painel_completo.png` — RMSD do backbone deve estabilizar
   (< 0,3 nm de variação no platô) para considerar o complexo estável.
2. Se o ligante dissociar (RMSD do ligante > 1 nm sem retorno, como ocorreu
   com XP352-GORE4), documentar e não prosseguir com MM-GBSA como resultado
   válido — o valor de ΔG ainda será calculado, mas não é interpretável.
3. Comparar `ΔG bind (TOTAL)` de DN773-GORE12T com os resultados equivalentes
   das outras isoformas (ACR157/QCL936/XP273) para ranquear afinidade.
4. Atualizar a memória do projeto (`project_md_gromacs.md`) com o resultado.
