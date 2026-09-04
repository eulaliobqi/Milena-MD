# Interação e afinidade de ligação da tripsina DN2954 com os peptídeos GORE1-2T e GORE1-2T(GGS)3: preparação, protonação em pH 8,2 e geração de pose por co-folding (Boltz-2)

*Bloco 1 de `novas-analises` — etapa pré-dinâmica molecular. Numeração de
resíduo do receptor sempre na convenção original do PDB (35–266;
His79/Asp126/Ser222/Asp216), a mesma já usada em
`assets/samplesheet_dn2954_gore12t.csv`. Escopo corrigido em 2026-09-03:
este documento mostra apenas o que segue para a próxima etapa (dinâmica
molecular em triplicata) — a pose gerada por co-folding com **Boltz-2**.
O docking guiado por restrições com **HADDOCK3** continua descrito nos
Métodos e usado como validação cruzada independente do sítio de ligação
(ele concorda com a região de interface e a tétrade catalítica), mas sua
própria pose de peptídeo não é usada em nenhuma etapa seguinte — o
motivo (limite de escopo do refinamento `flexref`, não um erro pontual)
está resumido na seção 2.2 e investigado por completo em
`comparativo_vs_DN2954-GORE12T.md` (seção 7), documento que preserva a
comparação completa entre os dois métodos para fins de proveniência.*

## 1. Material e Métodos

### 1.1 Sistema e dados de entrada

Receptor: tripsina digestiva **DN2954** (232 resíduos, His79/Asp126/
Ser222/Asp216 como tétrade catalítica clássica intacta), já modelada e
usada em trabalho anterior (`DN2954-GORE12T/`). Dois ligantes peptídicos:
**GORE1-2T** (21 aa, repetição `VLRVLKVLRVLKVLRVLKVLR`, sem estrutura
prévia) e **GORE1-2T(GGS)3** (75 aa, mesma repetição com três
espaçadores `(GGS)3` entre unidades) — sequência idêntica, resíduo a
resíduo, ao ligante `GORE12T` já usado e totalmente caracterizado por
100 ns de dinâmica molecular em `DN2954-GORE12T/`, o que permite usá-lo
como comparação direta de consistência de método (seção 2.7).

A identidade de sequência do receptor reaproveitado (`data/DN2954-
receptor.pdb`) contra `DN2954.fa` foi confirmada por extração da
sequência Cα (Biopython 1.83; Cock et al. 2009); faltam apenas os 3
resíduos C-terminais (`TIV`), região desordenada não resolvida no
modelo estrutural original.

### 1.2 Ensemble conformacional do peptídeo (entrada para o HADDOCK3)

Para cada peptídeo, três conformeros iniciais foram construídos a
partir da sequência com **PyMOL 3.1.0** (Schrödinger, LLC — The PyMOL
Molecular Graphics System), comando `cmd.fab()`: hélice-α (`ss=1`),
fita-β antiparalela estendida (`ss=2`) e poliprolina-II (construída
como cadeia estendida `ss=0` seguida de fixação por resíduo dos
ângulos diedrais canônicos φ=−75°, ψ=145° via `cmd.set_dihedral`). Os
três modelos foram combinados num único PDB multi-modelo com
`pdb_mkensemble` (pdb-tools 2.7.0; Rodrigues et al. 2018) e validados
com `pdb_chkensemble`. Este ensemble multi-confôrmero, em vez de uma
pose única de peptídeo, segue a recomendação corrente do próprio
HADDOCK3 para docking de peptídeos flexíveis (Giulini et al. 2025) —
usado como entrada do docking descrito em 1.3.

### 1.3 Validação cruzada do sítio de ligação (HADDOCK3)

Docking guiado por restrições ambíguas de interação (AIR) com
**HADDOCK3 versão 2026.7.0** (Giulini et al. 2025), executado
localmente (`mode = "local"`, `ncores = 12`). Usado neste bloco como
**confirmação independente de onde o peptídeo liga** (região de
interface, tétrade catalítica) — não como fonte da pose de MD; motivo
na seção 2.2.

- **Restrições:** resíduo ativo do receptor = tétrade catalítica
  (His79, Asp126, Ser222, Asp216); shell passivo calculado por
  `haddock3-restraints passive_from_active` (raio padrão 6,5 Å),
  retornando 31 resíduos de superfície ao redor do sítio ativo. Peptídeo
  tratado como totalmente passivo (todos os resíduos, sem restrição
  ativa), refletindo ausência de hipótese prévia sobre o resíduo P1 de
  contato. Tabela ambígua gerada por
  `haddock3-restraints active_passive_to_ambig`.
- **Fluxo (`.toml`)**, adaptado de uma corrida HADDOCK3 anterior do
  grupo já validada (`samera-wwp1-mbd2/03_docking/run_wwp1-mbd2.toml`):
  `topoaa` → `rigidbody` (`sampling = 1000`, `ambig_fname` = tabela AIR)
  → `caprieval` → `seletop` (`select = 200`) → `flexref` (peptídeo
  inteiro declarado semi-flexível, `nseg2=1`; ver seção 2.2) →
  `caprieval` → `emref` → `caprieval` → `clustfcc` (parâmetros padrão) →
  `seletopclusts` → `caprieval` final.
- Pose de referência de cada sistema: melhor modelo do cluster de
  menor (mais negativo) `HADDOCK score` no `caprieval` final — usada
  apenas para os números de validação citados na seção 2.2 (cluster
  score, recuperação de interface), não como estrutura de partida de MD.

### 1.4 Geração de pose — co-folding (Boltz-2, método adotado)

Co-folding direto receptor+peptídeo a partir das sequências, sem
posicionamento nem restrição prévia, com **Boltz versão 2.2.1**
(Passaro et al. 2025), ambiente `boltz2-env` dedicado. **Esta é a pose
usada como estrutura inicial da próxima etapa (MD em triplicata)** —
`model_0` de cada sistema.

- Entrada: YAML `version: 1`, duas entradas `protein` (cadeias A e B).
- Parâmetros: `boltz predict --use_msa_server --no_kernels
  --diffusion_samples 3 --output_format pdb`. MSA obtido do servidor
  público ColabFold (`api.colabfold.com`), pareamento `greedy`.
  `--no_kernels` usado por compatibilidade de kernel CUDA com a GPU
  RTX 5070 Ti (Blackwell); não afeta a arquitetura do modelo.
- **Não solicitado o bloco `properties: affinity`**: esse módulo do
  Boltz-2 só aceita ligante de molécula pequena (`ligand`/SMILES) como
  `binder`, não uma cadeia `protein`/peptídica — confirmado
  empiricamente (falha `ValueError: Chain B is not a ligand!`) antes de
  descartar o recurso, não assumido de antemão. A avaliação de
  confiança de cada amostra usou portanto só as métricas nativas de
  estrutura do Boltz-2 (`confidence_score`, `iptm`, `ptm`,
  `complex_plddt`, de `confidence_*.json`, e o pLDDT por resíduo
  gravado no campo B-factor do PDB de saída).
- Três amostras de difusão por sistema (não uma única), seguindo a
  recomendação de checagem de consistência de pose entre amostras já
  registrada em `docs/03_metodologia_padrao_ouro.md` (motivada por
  Kim et al. 2026 — confiança alta do modelo não garante pose
  fisicamente correta). Essa checagem de consistência foi decisiva:
  ver seção 2.2/2.7.

### 1.5 Controle de qualidade estrutural

Etapa nova em relação ao pipeline de produção do Milena-MD (que hoje
não tem controle de qualidade estrutural formal). Sobre cada uma das 3
amostras de difusão Boltz-2 de cada sistema, a geometria da díade
catalítica foi calculada estaticamente com um script próprio
(`scripts/triad_distances.py`, Biopython 1.83 `Bio.PDB`), medindo as
distâncias heavy-atom His-NE2···Ser-OG e His-ND1···Asp-OD1/OD2, com
limiar de aceitação ≤4,5 Å (geometria de par catalítico plausível). A
mesma checagem foi aplicada à pose HADDOCK (validação cruzada, seção
1.3); os números completos das duas ferramentas lado a lado estão em
`comparativo_vs_DN2954-GORE12T.md`.

### 1.6 Protonação em pH 8,2 e avaliação de resíduos ionizáveis

Protonação com **PDB2PQR versão 3.7.1** (Dolinsky et al. 2004) usando
o motor de estados de titulação **PROPKA versão 3.5.1** (Olsson et al.
2011; Søndergaard et al. 2011): `pdb2pqr30 --ff AMBER
--titration-state-method propka --with-ph 8.2 --keep-chain
--pdb-output`, pH fixado no valor já padronizado pelo grupo para
intestino médio de lepidóptero. A tabela de pKa por resíduo, impressa
pelo PROPKA em STDERR (não STDOUT), foi extraída por script próprio
(`scripts/protonate_ph82.py`), filtrando especificamente HIS, ASP, GLU
e CYS e sinalizando qualquer resíduo cujo pKa previsto cruzasse o pH
de referência (8,2) na direção incomum (HIS>8,2 → carregado; ASP/GLU>
8,2 → protonado; CYS com pKa próximo de 8,2 → estado incerto). Aplicado
à pose Boltz-2 `model_0` de cada sistema (resultados na seção 2.4); a
mesma checagem na pose HADDOCK está em `comparativo_vs_DN2954-GORE12T.md`
seção 4.

### 1.7 Análise de interação (interface receptor-peptídeo)

Para cada pose protonada, resíduos do receptor (cadeia A) com pelo
menos um átomo pesado a ≤4,5 Å de qualquer átomo pesado do peptídeo
(cadeia B) foram listados com `NeighborSearch` do Biopython
(`scripts/interface_contacts.py`) — equivalente de pose única, sem
trajetória de MD, ao mapa de contato por frequência já usado em
produção (`bin/contact_map.py`).

### 1.8 Estimativa de afinidade de ligação

Energia livre de ligação estimada com **PRODIGY-prot versão 2.3.0**
(Xue et al. 2016), método baseado em contatos de interface calibrado
em complexos proteína-proteína estruturados: `prodigy <pdb>
--selection A B`, temperatura padrão 25 °C, sobre a pose Boltz-2
`model_0` já protonada em pH 8,2. Ressalva declarada: a validade
estatística do PRODIGY para uma cadeia B curta/não-globular (peptídeo,
não domínio proteico) não está comprovada — reportado do mesmo modo
que em uso anterior do grupo para o peptídeo GORE3 (5 aa) em
`eulalio-pos-doc/codigo/fase9_blocoF6c_docking_energies/run_prodigy.py`,
com a mesma ressalva propagada.

---

## 2. Resultados

### 2.1 A estrutura do receptor e do ligante GGS3 são reaproveitáveis sem retrabalho

A sequência Cα extraída de `DN2954-receptor.pdb` reproduz exatamente
`DN2954.fa` (232/232 resíduos, faltando apenas os 3 resíduos
C-terminais desordenados `TIV`), e a sequência de
`GORE12T-ligand-DN2954.pdb` é idêntica, resíduo a resíduo, a
`GORE1-2T(GGS)3.fa` (75/75). Nenhuma nova predição de estrutura de
receptor foi necessária; o ligante GGS3 pôde ser tratado como o mesmo
objeto molecular já validado por MD, permitindo comparação direta de
método em vez de comparação entre moléculas diferentes.

### 2.2 Boltz-2 converge em pose única e bem definida para os dois peptídeos; HADDOCK3 confirma o sítio mas não fornece uma pose de peptídeo independente

O co-folding Boltz-2 concorda em confiança alta para as 3 amostras de
difusão de ambos os sistemas (`confidence_score` 0,86–0,90; Fig. 1),
com uma diferença sistemática entre peptídeos: **o peptídeo curto
GORE1-2T é predito com confiança de interface (ipTM) mais alta e mais
estável entre amostras (0,81/0,83/0,67) do que o construto ligado
GORE1-2T(GGS)3 (0,73/0,58/0,46)** — consistente com a maior
flexibilidade conformacional esperada de um peptídeo 75 aa com três
espaçadores `(GGS)3` frente a um peptídeo curto de 21 aa.

**Validação cruzada:** o docking guiado por restrições HADDOCK3 (seção
1.3) concorda com essa região de ligação — cluster líder claramente
separado dos demais por `HADDOCK score` nos dois sistemas (GORE1-2T:
−67,87 ± 3,53 vs. segundo colocado −55,36 ± 4,52; GORE1-2T(GGS)3:
−56,32 ± 1,89 vs. segundo colocado −53,33 ± 3,04) e recuperação de
interface comparável à do Boltz-2 (seção 2.5). Essa concordância é
usada aqui só como confirmação de que a região de ligação identificada
é robusta a dois métodos independentes, **não** como fonte de uma pose
alternativa de peptídeo: investigação dirigida (RMSD de backbone,
detalhe completo em `comparativo_vs_DN2954-GORE12T.md` seção 7) mostrou
que o backbone do peptídeo na pose HADDOCK tem RMSD de apenas
0,38–0,43 Å contra o conformero hélice-α de entrada (seção 1.2) — ou
seja, o `flexref` do HADDOCK3 não o alterou de fato. Mesmo depois de
reconfigurar `flexref` para tratar o peptídeo inteiro como
semi-flexível e reprocessar a partir dessa etapa (`--restart 4`,
reaproveitando o rigidbody já calculado), o `HADDOCK score` não mudou
(estatisticamente idêntico) mas a conformação também não: RMSD
0,41–0,75 Å contra a mesma hélice, só reposicionada em relação ao
receptor. Conclusão: **não é um parâmetro mal ajustado, é um limite de
escopo do `flexref`** — refinamento local por SA de torção em rajadas
curtas, sem capacidade de desfazer uma α-hélice já formada. Por isso a
pose usada na próxima etapa é sempre a do **Boltz-2 `model_0`** (Fig.
6) — HADDOCK3 fica como confirmação independente do sítio, não como
gerador de pose.

### 2.3 A tétrade catalítica permanece geometricamente intacta nas 6 poses Boltz-2 avaliadas

As distâncias His-NE2···Ser-OG (3,19–3,30 Å) e His-ND1···Asp-OD1/OD2
(2,75–2,83 Å, menor das duas) ficaram dentro da faixa de par catalítico
plausível em **todas** as 3 amostras de difusão dos dois peptídeos
(Fig. 2). Nenhuma amostra foi descartada por tríade rompida. A
distância Ser-OG···Asp216 (9,7–10,7 Å em todas as poses) não constitui
achado de baixa qualidade: Asp216 é o resíduo de especificidade do
fundo do bolsão S1, que interage com a cadeia lateral do resíduo P1 do
peptídeo ligante, não faz ponte de hidrogênio direta com a serina
catalítica.

### 2.4 A protonação em pH 8,2 confirma o padrão esperado de carga, sem anomalias na pose adotada

Nas 2 poses protonadas (Boltz-2 `model_0`, os dois peptídeos), todo ASP
e GLU (17 resíduos por sistema) foi previsto carregado (ASP 2,89–5,03;
GLU 2,63–5,32; ambos bem abaixo de 8,2; Fig. 3), e todas as 6 cisteínas
por sistema retornaram pKa=99,99 do PROPKA — sentinela de resíduo
não-titulável, consistente com a arquitetura de pontes dissulfeto já
conhecida da tripsina, sem indício de tiolato livre em pH 8,2. Todas as
6 histidinas por sistema, incluindo a His79 catalítica (pKa 6,97 no
GORE1-2T; 7,04 no GGS3), saem neutras (HID/HIE, pKa 4,15–7,04) — **sem
exceção** nas poses adotadas. (A pose HADDOCK descartada mostrava uma
His79 anômala, pKa 8,55, especificamente no GORE1-2T — artefato de
ponte salina com o C-terminal livre do peptídeo naquela pose
específica, sem relação com a pose Boltz-2 usada aqui; causa raiz
completa em `comparativo_vs_DN2954-GORE12T.md` seção 4.)

### 2.5 A interface de ligação prevista recupera a maior parte dos contatos já validados por 100 ns de dinâmica molecular

Os resíduos de interface (≤4,5 Å) previstos pelo Boltz-2 para o
GORE1-2T(GGS)3 — molécula idêntica à já simulada em `DN2954-GORE12T/`
— recuperam 9 dos ~20 resíduos de maior frequência de contato
observados ao longo de 100 ns de MD (55, 56, 63, 120, 218, 219, 220,
240, 242, entre outros; Fig. 5), apesar de o método não ter usado essa
trajetória como referência ou restrição. O peptídeo novo GORE1-2T
(sem espaçadores, sem histórico de MD) recupera um subconjunto
igualmente numeroso (9 resíduos em comum), concentrado na mesma alça
219–222/237–242 ao redor do sítio ativo — o mesmo bolsão de ligação,
não um sítio alternativo. O docking HADDOCK3 (validação cruzada,
seção 2.2) recupera um conjunto parcialmente sobreposto (4–9 resíduos
em comum, números completos em `comparativo_vs_DN2954-GORE12T.md`
seção 5), consistente com a mesma região.

### 2.6 A afinidade prevista é favorável para os dois peptídeos

PRODIGY-prot prevê ligação favorável para a pose adotada nos dois
sistemas: GORE1-2T ΔG = −10,6 kcal/mol (Kd ≈ 1,5×10⁻⁸ M); GORE1-2T(GGS)3
ΔG = −11,2 kcal/mol (Kd ≈ 5,9×10⁻⁹ M) — Fig. 4, diferença de ~0,6
kcal/mol entre os dois peptídeos, não decisiva isoladamente. O valor de
energia livre já existente para `DN2954-GORE12T` (estimativa de
Interaction Entropy pós-MD, −104,65 kcal/mol) não é comparável em
magnitude absoluta a este valor de PRODIGY — é uma estimativa em
vácuo, sem termo de solvatação, que o próprio grupo já classifica como
instável (desvio padrão >10 kT); a comparação válida entre as duas
gerações de análise é de direção (ligação favorável em ambas), não de
módulo.

### 2.7 Síntese: pose adotada e ressalva remanescente para a MD

O peptídeo GORE1-2T(GGS)3 — reprocessado do zero com um protocolo que
não teve acesso à trajetória de MD anterior — converge com o resultado
já validado de `DN2954-GORE12T` em geometria de sítio ativo, região de
interface e sinal de afinidade favorável, usando a pose gerada por
co-folding (Boltz-2) e com validação independente adicional do
HADDOCK3 no próprio sítio de ligação. O peptídeo novo GORE1-2T, sem
histórico prévio, liga na mesma região do sítio ativo com afinidade
favorável.

**Estrutura adotada para a próxima etapa (dinâmica molecular em
triplicata): a pose Boltz-2 `model_0`, para os dois sistemas** (Fig.
6). **Ressalva registrada, não escondida:** para o GORE1-2T
especificamente, as 3 amostras de difusão do próprio Boltz-2 divergem
entre si na conformação livre do peptídeo (distância ponta-a-ponta
22,0–36,6 Å; pLDDT médio 47–56, abaixo do limiar costumeiro de
confiança) — ao contrário do GGS3, cujas 3 amostras convergem
fortemente (16,0–16,6 Å; pLDDT 62–75). Isso não é um problema da região
de interface (que recupera contatos validados nos dois peptídeos,
seção 2.5) — é uma incerteza real sobre a conformação da porção livre
do GORE1-2T curto, a esclarecer pela própria MD de produção (por
exemplo comparando as réplicas da triplicata), não algo já resolvido
nesta etapa de docking/co-folding. Investigação completa desta seção,
incluindo a pose HADDOCK descartada e a tentativa (sem sucesso) de
corrigi-la reconfigurando o `flexref`: `comparativo_vs_DN2954-GORE12T.md`.

---

## Legendas das figuras

Todas as figuras: PDF vetorial (fontes TrueType embutidas) + PNG a 600
dpi, estilo de partida "Nature" (`scientific-visualization`/
`high-impact-figures`, snapshot de 2026-07-23 — ponto de partida
visual, não uma certificação de conformidade de submissão), paleta
Okabe-Ito sobre fundo branco. Em toda figura a categoria (sistema,
método, resíduo) já está identificada por posição/rótulo do eixo, não
só por cor. Fontes de dados e manifestos de exportação
(`.export.json`, com proveniência) em `novas-analises/figuras/`. Todas
as figuras abaixo usam exclusivamente a pose Boltz-2 (a que segue para
a MD); a comparação completa com a pose HADDOCK descartada está
arquivada em `figuras/arquivo_divergencia_haddock_boltz/` (fora da
numeração deste artigo) e em `comparativo_vs_DN2954-GORE12T.md`.

**Fig. 1 — Boltz-2: métricas de confiança por amostra de difusão.**
`confidence_score`, ipTM, pTM e pLDDT do complexo, lidos diretamente de
`confidence_*.json`, para as 3 amostras (`--diffusion_samples 3`) de
cada sistema.

**Fig. 2 — Geometria da díade catalítica His79/Ser222/Asp126, poses
Boltz-2.** Distâncias heavy-atom His-NE2···Ser-OG e His-ND1···Asp-OD
(menor das duas, OD1/OD2) em cada uma das 3 amostras de difusão × 2
peptídeos (6 poses). Linha tracejada = limiar de aceitação (4,5 Å)
usado no controle de qualidade.

**Fig. 3 — pKa previsto de HIS/ASP/GLU nas 2 poses Boltz-2 `model_0`
protonadas em pH 8,2.** PROPKA 3.5.1 via PDB2PQR 3.7.1. Linha tracejada
= pH de referência (8,2). Nenhuma anomalia sinalizada nas poses
adotadas. CYS omitida do painel — pKa=99,99 (sentinela de
não-titulável/dissulfeto) em todas as 6 cisteínas das 2 poses.

**Fig. 4 — Afinidade prevista (PRODIGY-prot 2.3.0), pose Boltz-2
`model_0` protonada em pH 8,2.** ΔG por `prodigy <pdb> --selection A B`
(25 °C), para os 2 peptídeos.

**Fig. 5 — Recuperação da interface já validada por 100 ns de MD, pose
Boltz-2.** Nº de resíduos do receptor (contato ≤4,5 Å) em comum com os
20 resíduos de maior frequência de contato na MD de 100 ns de
`DN2954-GORE12T` (`analise_extra/interface_residues.csv`,
`max_contact_freq`), por sistema. Linha tracejada = tamanho total do
conjunto de referência (n=20).

**Fig. 6 — Pose Boltz-2 `model_0` adotada como estrutura inicial de
MD, renderização estrutural.** PyMOL 3.1.0, cadeia A (receptor) de
cada pose superposta rigidamente (`cmd.align`) sobre
`data/DN2954-receptor.pdb`; tétrade catalítica em verde (sticks).
**(a)** DN2954 × GORE1-2T (21 aa) e **(b)** DN2954 ×
GORE1-2T(GGS)3 (75 aa): mesma câmera/escala nos dois painéis. Script de
renderização: `scripts/render_boltz_poses.py` (executado no servidor,
ambiente `structure`); composição: `scripts/make_fig6_boltz_pose.py`.

## Referências

Cock, P. J. A. et al. Biopython: freely available Python tools for
computational molecular biology and bioinformatics. *Bioinformatics*
**25**, 1422–1423 (2009).

Schrödinger, LLC. *The PyMOL Molecular Graphics System*, Version 3.1
(2025).

Rodrigues, J. P. G. L. M., Teixeira, J. M. C., Trellet, M. & Bonvin, A.
M. J. J. pdb-tools: a swiss army knife for molecular structures.
*F1000Research* **7**, 1961 (2018).

Giulini, M. et al. HADDOCK3. *J. Chem. Inf. Model.* (2025).
doi:10.1021/acs.jcim.5c00969

Passaro, S. et al. Boltz-2: Towards Accurate and Efficient Binding
Affinity Prediction. *bioRxiv* (2025). PMC12262699.

Kim, J. et al. Large-scale evaluation of co-folding pose accuracy vs.
model confidence. (2026). [citado em
`docs/03_metodologia_padrao_ouro.md`, referência interna `kim2026largescale`]

Dolinsky, T. J., Nielsen, J. E., McCammon, J. A. & Baker, N. A.
PDB2PQR: an automated pipeline for the setup of Poisson-Boltzmann
electrostatics calculations. *Nucleic Acids Res.* **32**, W665–W667
(2004).

Olsson, M. H. M., Søndergaard, C. R., Rostkowski, M. & Jensen, J. H.
PROPKA3: Consistent Treatment of Internal and Surface Residues in
Empirical pKa Predictions. *J. Chem. Theory Comput.* **7**, 525–537
(2011).

Søndergaard, C. R., Olsson, M. H. M., Rostkowski, M. & Jensen, J. H.
Improved Treatment of Ligands and Coupling Effects in Empirical
Calculation and Rationalization of pKa Values. *J. Chem. Theory
Comput.* **7**, 2284–2295 (2011).

Xue, L. C., Rodrigues, J. P., Kastritis, P. L., Bonvin, A. M. &
Vangone, A. PRODIGY: a web server for predicting the binding affinity
of protein-protein complexes. *Bioinformatics* **32**, 3676–3678
(2016).

## Reprodutibilidade — localização de código e dados

Scripts: `novas-analises/scripts/{protonate_ph82.py, triad_distances.py,
run_prodigy_batch.py, interface_contacts.py}`,
`novas-analises/build_peptide_ensemble.py`. Configurações HADDOCK3
(validação cruzada, seção 1.3):
`novas-analises/DN2954-GORE1-2T{,-GGS3}/haddock/run_dn2954-gore1-2t{,-ggs3}.toml`
(`[flexref] nseg2=1` — configuração corrigida; run original preservado
em `haddock/run1-guided_orig_frozen-peptide/`). Entradas Boltz-2 (pose
adotada): `novas-analises/DN2954-GORE1-2T{,-GGS3}/boltz2/*.yaml`.
Resultados brutos e tabelas por sistema:
`novas-analises/DN2954-GORE1-2T{,-GGS3}/{qualidade,protonacao_ph8.2,interacao,afinidade}/`
(contêm ambos os métodos; as figuras deste artigo filtram só a linha/
pose Boltz-2 — ver `scripts/make_figures_pub.py`). Comparação completa
com a MD anterior e com a pose HADDOCK descartada:
`novas-analises/comparativo_vs_DN2954-GORE12T.md`. Figuras 1–5
(PDF+PNG 600dpi, com manifesto `.export.json`) e script gerador:
`novas-analises/figuras/`, `novas-analises/scripts/make_figures_pub.py`.
Fig. 6 (pose Boltz-2, renderização estrutural): script PyMOL
`novas-analises/scripts/render_boltz_poses.py` (roda no servidor,
ambiente `structure`, produz `figuras/_render/boltz_pose_*.png`) +
script de composição `novas-analises/scripts/make_fig6_boltz_pose.py`
(roda localmente, aplica autocrop + estilo `scientific-visualization`).
Divergência de conformação HADDOCK vs. Boltz-2 (seção 2.2/2.7, detalhe
em `comparativo_vs_DN2954-GORE12T.md` seção 7): números reproduzidos
por `novas-analises/scripts/diagnose_pose_divergence.py`, sobre os PDBs
em `novas-analises/diagnostico_divergencia_pose/`; render/figura
arquivados em `figuras/arquivo_divergencia_haddock_boltz/` (fora da
numeração deste artigo).
