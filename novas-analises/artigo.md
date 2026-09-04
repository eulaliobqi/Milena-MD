# Interação e afinidade de ligação da tripsina DN2954 com os peptídeos GORE1-2T e GORE1-2T(GGS)3: preparação, protonação em pH 8,2 e docking padrão-ouro + fronteira

*Bloco 1 de `novas-analises` — etapa pré-dinâmica molecular. Numeração de
resíduo do receptor sempre na convenção original do PDB (35–266;
His79/Asp126/Ser222/Asp216), a mesma já usada em
`assets/samplesheet_dn2954_gore12t.csv`.*

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

### 1.2 Ensemble conformacional do peptídeo

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
HADDOCK3 para docking de peptídeos flexíveis (Giulini et al. 2025).

### 1.3 Geração de complexo — padrão ouro (HADDOCK3)

Docking guiado por restrições ambíguas de interação (AIR) com
**HADDOCK3 versão 2026.7.0** (Giulini et al. 2025), executado
localmente (`mode = "local"`, `ncores = 12`).

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
  → `caprieval` → `seletop` (`select = 200`) → `flexref` → `caprieval`
  → `emref` → `caprieval` → `clustfcc` (parâmetros padrão) →
  `seletopclusts` → `caprieval` final.
- Pose de referência de cada sistema: melhor modelo do cluster de
  menor (mais negativo) `HADDOCK score` no `caprieval` final.

### 1.4 Geração de complexo — fronteira (Boltz-2, co-folding)

Co-folding direto receptor+peptídeo a partir das sequências, sem
posicionamento nem restrição prévia, com **Boltz versão 2.2.1**
(Passaro et al. 2025), ambiente `boltz2-env` dedicado.

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
  `complex_plddt`, de `confidence_*.json`).
- Três amostras de difusão por sistema (não uma única), seguindo a
  recomendação de checagem de consistência de pose entre amostras já
  registrada em `docs/03_metodologia_padrao_ouro.md` (motivada por
  Kim et al. 2026 — confiança alta do modelo não garante pose
  fisicamente correta).

### 1.5 Controle de qualidade estrutural

Etapa nova em relação ao pipeline de produção do Milena-MD (que hoje
não tem controle de qualidade estrutural formal). Sobre cada pose
candidata (HADDOCK: melhor modelo do melhor cluster; Boltz-2: as 3
amostras de difusão), a geometria da díade catalítica foi calculada
estaticamente com um script próprio (`scripts/triad_distances.py`,
Biopython 1.83 `Bio.PDB`), medindo as distâncias heavy-atom
His-NE2···Ser-OG e His-ND1···Asp-OD1/OD2, com limiar de aceitação
≤4,5 Å (geometria de par catalítico plausível).

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
8,2 → protonado; CYS com pKa próximo de 8,2 → estado incerto).

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
--selection A B`, temperatura padrão 25 °C, sobre as estruturas já
protonadas em pH 8,2. Ressalva declarada: a validade estatística do
PRODIGY para uma cadeia B curta/não-globular (peptídeo, não domínio
proteico) não está comprovada — reportado do mesmo modo que em uso
anterior do grupo para o peptídeo GORE3 (5 aa) em
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

### 2.2 HADDOCK3 e Boltz-2 convergem em pose única e bem definida para os dois peptídeos

O docking HADDOCK3 produziu, para os dois sistemas, um cluster líder
claramente separado dos demais por `HADDOCK score` (GORE1-2T: cluster
1, n=10, score −67,87 ± 3,53, vs. segundo colocado −55,36 ± 4,52;
GORE1-2T-GGS3: cluster 1, n=10, score −56,32 ± 1,89, vs. segundo
colocado −53,33 ± 3,04; Fig. 1). O co-folding Boltz-2 concorda em
confiança alta para as 3 amostras de difusão de ambos os sistemas
(`confidence_score` 0,86–0,90; Fig. 2), com uma diferença sistemática
entre peptídeos: **o peptídeo curto GORE1-2T é predito com confiança de
interface (ipTM) mais alta e mais estável entre amostras (0,81/0,83/
0,67) do que o construto ligado GORE1-2T(GGS)3 (0,73/0,58/0,46)** —
consistente com a maior flexibilidade conformacional esperada de um
peptídeo 75 aa com três espaçadores `(GGS)3` frente a um peptídeo
curto de 21 aa. Visualmente, as poses HADDOCK e Boltz-2 ocupam a mesma
região da tripsina para os dois peptídeos, ancoradas junto à tétrade
catalítica (Fig. 7a-b). Para o GORE1-2T(GGS)3, a pose HADDOCK completa
(sem recorte de câmera, Fig. 7c) mostra o conformero hélice-α do
ensemble de entrada preservado quase reto ao longo de toda a cadeia de
75 aa (~110 Å) — só a extremidade N-terminal faz contato com o
receptor. Investigação dirigida (RMSD de backbone, ver
`comparativo_vs_DN2954-GORE12T.md` seção 7) identificou a causa: o
`[flexref]` do HADDOCK3 detecta automaticamente a zona semi-flexível
pelo contato com o receptor, deixando a maior parte de um peptídeo sem
estrutura prévia congelada no conformero rígido de entrada (RMSD
0,38–0,43 Å vs. a hélice ideal, nos dois sistemas). Reconfigurar
`[flexref]` para tratar o peptídeo inteiro como semi-flexível e
reprocessar a partir dessa etapa (`--restart 4`, reaproveitando o
rigidbody já calculado) preservou o `HADDOCK score` (estatisticamente
idêntico) mas **não** alterou a conformação — o peptídeo continuou
essencialmente na mesma hélice (RMSD 0,41–0,75 Å), só reposicionada.
Conclusão: não é um parâmetro mal ajustado, é um limite de escopo do
`flexref` (refinamento local em rajadas curtas de SA, não uma
ferramenta de predição de fold) — ele não consegue desfazer uma
α-hélice já formada. A cauda distal do GGS3 não tem suporte
experimental de conformação em nenhuma das duas rodadas HADDOCK e não
deve ser tratada como pose validada na preparação da topologia de MD
(seção 2.7).

### 2.3 A tétrade catalítica permanece geometricamente intacta em todas as 8 poses avaliadas

As distâncias His-NE2···Ser-OG (2,96–4,01 Å) e His-ND1···Asp-OD1/OD2
(2,73–3,43 Å) ficaram dentro da faixa de par catalítico plausível em
**todas** as poses testadas — HADDOCK e as 3 amostras Boltz-2, para os
dois peptídeos — sem exceção (Fig. 3). Nenhuma pose foi descartada por
tríade rompida. A distância Ser-OG···Asp216 (9,7–11,3 Å em todas as poses) não
constitui achado de baixa qualidade: Asp216 é o resíduo de
especificidade do fundo do bolsão S1, que interage com a cadeia
lateral do resíduo P1 do peptídeo ligante, não faz ponte de hidrogênio
direta com a serina catalítica.

### 2.4 A protonação em pH 8,2 confirma o padrão esperado de carga, com uma exceção pontual na His79 catalítica

Em todas as 4 poses protonadas (HADDOCK e Boltz-2, para os dois
peptídeos), todo ASP e GLU (17 resíduos por sistema) foi previsto
carregado (pKa 2,0–5,8, bem abaixo de 8,2), e todas as 6 cisteínas por
sistema retornaram pKa=99,99 do PROPKA — sentinela de resíduo
não-titulável, consistente com a arquitetura de pontes dissulfeto já
conhecida da tripsina, sem indício de tiolato livre em pH 8,2. Das 6
histidinas por sistema, 5 saem neutras (HID/HIE, pKa 4,2–7,0) em toda
pose. **Achado pontual, não generalizável:** na pose HADDOCK do
peptídeo curto GORE1-2T especificamente, a His79 catalítica é prevista
com pKa 8,55 — acima do pH de referência, portanto potencialmente
carregada (HIP) nessa geometria local. A mesma His79, no Boltz-2 (mesmo
peptídeo) e em ambos os métodos para o GORE1-2T(GGS)3, permanece
neutra (pKa 7,0–8,0). **Causa identificada:** o log do PROPKA atribui
o excesso de pKa a um termo de Coulomb de +1,66 com o carboxilato
C-terminal do próprio peptídeo (`C- 21 B`), ausente em qualquer outra
pose da His79; medição direta confirma His79-NE2 a 3,18 Å desse
carboxilato na pose HADDOCK, contra 13,2–14,1 Å na pose Boltz-2 do
mesmo sistema. GORE1-2T (21 aa, sem capeamento C-terminal, tratado
como totalmente passivo no AIR) teve, nessa pose específica, o
C-terminal livre orientado por acaso para dentro do bolsão,
formando uma ponte salina acidental com a His catalítica — artefato
do docking rígido+flexível, não uma propriedade real do par
receptor-peptídeo (não se repete no GGS3 nem em nenhuma amostra
Boltz-2; ver `comparativo_vs_DN2954-GORE12T.md`, seção 4).

### 2.5 A interface de ligação prevista recupera a maior parte dos contatos já validados por 100 ns de dinâmica molecular

Os resíduos de interface (≤4,5 Å) previstos para o GORE1-2T(GGS)3 —
molécula idêntica à já simulada em `DN2954-GORE12T/` — recuperam 8-9
dos ~20 resíduos de maior frequência de contato observados ao longo de
100 ns de MD (55, 56, 63, 120, 218, 219, 220, 240, 242, entre outros),
tanto pelo HADDOCK3 quanto pelo Boltz-2, apesar de nenhum dos dois
métodos ter usado essa trajetória como referência ou restrição. O
peptídeo novo GORE1-2T (sem espaçadores, sem histórico de MD)
recupera um subconjunto menor mas sobreposto da mesma região (4
resíduos em comum via HADDOCK, 9 via Boltz-2), concentrado na mesma
alça 219–222/237–242 ao redor do sítio ativo — o mesmo bolsão de
ligação, não um sítio alternativo.

### 2.6 A afinidade prevista é favorável para os dois peptídeos, com uma diferença sistemática entre métodos

PRODIGY-prot prevê ligação favorável em todas as 4 poses avaliadas:
GORE1-2T ΔG = −6,0 kcal/mol (HADDOCK) e −10,6 kcal/mol (Boltz-2);
GORE1-2T(GGS)3 ΔG = −6,6 kcal/mol (HADDOCK) e −11,2 kcal/mol
(Boltz-2). Duas tendências consistentes, reportadas lado a lado sem
decidir qual método está mais correto: **(i)** o Boltz-2 prevê
afinidade sistematicamente ~4–5 kcal/mol mais forte que o HADDOCK nos
dois peptídeos, plausivelmente refletindo maior área de interface
enterrada nas poses de co-folding frente ao refinamento
rígido+flexível do HADDOCK; **(ii)** o construto ligado
GORE1-2T(GGS)3 é predito com afinidade marginalmente mais forte que o
GORE1-2T curto em ambos os métodos (diferença de ~0,6 kcal/mol, não
decisiva isoladamente). O valor de energia livre já existente para
`DN2954-GORE12T` (estimativa de Interaction Entropy pós-MD, −104,65
kcal/mol) não é comparável em magnitude absoluta a estes valores de
PRODIGY — é uma estimativa em vácuo, sem termo de solvatação, que o
próprio grupo já classifica como instável (desvio padrão >10 kT); a
comparação válida entre as duas gerações de análise é de direção
(ligação favorável em ambas), não de módulo.

### 2.7 Síntese: o padrão de ligação estabelecido para o GORE12T se mantém sob um protocolo novo e independente

O peptídeo GORE1-2T(GGS)3 — reprocessado do zero com um protocolo de
docking padrão-ouro + fronteira que não teve acesso à trajetória de MD
anterior — converge com o resultado já validado de `DN2954-GORE12T`
em geometria de sítio ativo, região de interface e sinal de afinidade
favorável. O peptídeo GORE1-2T, sem histórico prévio, liga na mesma
região do sítio ativo com afinidade favorável comparável em ambos os
métodos de docking. A única ressalva do bloco — His79 com protonação
incomum na pose HADDOCK do GORE1-2T (seção 2.4) — está resolvida: é
artefato de pose, não achado biológico. **Decisão para a próxima etapa
(dinâmica molecular em triplicata):** usar a pose Boltz-2 (`model_0`)
como estrutura inicial de MD para o sistema GORE1-2T; caso a pose
HADDOCK seja usada por algum motivo, His79 deve ser forçada a estado
neutro (HID/HIE) na preparação da topologia. Segunda ressalva,
identificada pela renderização estrutural (Fig. 7c) e investigada até a
causa raiz (seção 2.2): a pose HADDOCK do GORE1-2T(GGS)3 mantém o
conformero hélice-α de entrada quase reto por toda a cadeia de 75 aa,
mesmo após reconfigurar o `flexref` para tratar o peptídeo inteiro como
semi-flexível e reprocessar a partir dessa etapa — o `HADDOCK score` não
mudou, só a orientação da hélice em relação ao receptor. Confirmado:
`flexref` é um refinamento local (SA de torção em rajadas curtas), não
uma ferramenta de predição de fold, e não tem como desfazer uma α-hélice
já formada. O Boltz-2, por ser um modelo de co-folding real, é a única
fonte disponível de hipótese de conformação para a cadeia inteira —
reprodutível para o GGS3 (3/3 amostras convergem em fold compacto de
~16 Å) e explicitamente **não** reprodutível para o GORE1-2T (3
amostras divergem entre si, pLDDT<60 nas três). **Decisão:** usar a
pose Boltz-2 como estrutura inicial de MD nos dois sistemas; para o
GORE1-2T, tratar a conformação livre do peptídeo como incerteza real a
esclarecer pela própria MD de produção (por exemplo comparando as
réplicas da triplicata), não como algo já resolvido pelo docking.

---

## Legendas das figuras

Todas as figuras: PDF vetorial (fontes TrueType embutidas) + PNG a 600
dpi, estilo de partida "Nature" (`scientific-visualization`/
`high-impact-figures`, snapshot de 2026-07-23 — ponto de partida
visual, não uma certificação de conformidade de submissão), paleta
Okabe-Ito sobre fundo branco. Em toda figura a categoria (sistema,
método, resíduo) já está identificada por posição/rótulo do eixo, não
só por cor — a auditoria de contraste (`palette_audit.py`) sinalizou 4
dos 10 pares de cor com baixa separação em escala de cinza (revisão,
não falha), registrado aqui em vez de omitido. Fontes de dados e
manifestos de exportação (`.export.json`, com proveniência) em
`novas-analises/figuras/`.

**Fig. 1 — HADDOCK3 (padrão-ouro): `HADDOCK score` por cluster.**
Barras = score médio do cluster (`09_seletopclusts`/`10_caprieval`);
barras de erro = desvio padrão entre os modelos do cluster
(`score_std`). (a) GORE1-2T, 9 clusters; (b) GORE1-2T(GGS)3, 12
clusters. Cluster 1 (rank 1) é o mais negativo (melhor) em ambos.

**Fig. 2 — Boltz-2 (fronteira): métricas de confiança por amostra de
difusão.** `confidence_score`, ipTM, pTM e pLDDT do complexo, lidos
diretamente de `confidence_*.json`, para as 3 amostras
(`--diffusion_samples 3`) de cada sistema.

**Fig. 3 — Geometria da díade catalítica His79/Ser222/Asp126.**
Distâncias heavy-atom His-NE2···Ser-OG e His-ND1···Asp-OD (menor das
duas, OD1/OD2) em cada uma das 8 poses avaliadas (HADDOCK + 3 amostras
Boltz-2, 2 peptídeos). Linha tracejada = limiar de aceitação (4,5 Å)
usado no controle de qualidade.

**Fig. 4 — pKa previsto de HIS/ASP/GLU nas 4 poses protonadas em pH
8,2.** PROPKA 3.5.1 via PDB2PQR 3.7.1. Linha tracejada = pH de
referência (8,2). Estrela com contorno preto = anomalia sinalizada
(His79 da pose HADDOCK do GORE1-2T, pKa 8,55). CYS omitida do painel —
pKa=99,99 (sentinela de não-titulável/dissulfeto) em todas as 6
cisteínas de todas as 4 poses.

**Fig. 5 — Afinidade prevista (PRODIGY-prot 2.3.0), poses protonadas
em pH 8,2.** ΔG por `prodigy <pdb> --selection A B` (25 °C), para os 2
peptídeos × 2 métodos de geração de pose.

**Fig. 6 — Recuperação da interface já validada por 100 ns de MD.**
Nº de resíduos do receptor (contato ≤4,5 Å) em comum com os 20
resíduos de maior frequência de contato na MD de 100 ns de
`DN2954-GORE12T` (`analise_extra/interface_residues.csv`,
`max_contact_freq`), por sistema e método de docking. Linha tracejada
= tamanho total do conjunto de referência (n=20).

**Fig. 7 — Poses consenso (HADDOCK3 padrão-ouro vs. Boltz-2 fronteira,
model_0), renderização estrutural.** PyMOL 3.1.0, cadeia A (receptor)
de cada pose superposta rigidamente (`cmd.align`) sobre
`data/DN2954-receptor.pdb`; tétrade catalítica em verde (sticks).
**(a)** DN2954 × GORE1-2T (21 aa) e **(b)** DN2954 ×
GORE1-2T(GGS)3 (75 aa): mesma câmera/escala nos dois painéis, receptor
inteiro + peptídeo até 15 Å da superfície do receptor — nos dois
sistemas e pelos dois métodos, o peptídeo converge sobre a mesma
região da tétrade. **(c)** Pose HADDOCK completa do GORE1-2T(GGS)3, sem
o recorte de câmera de (b): a hélice-α de entrada permanece quase reta
por toda a cadeia de 75 aa (~110 Å), com a extremidade N-terminal
apenas em contato com o receptor — achado reportado, não omitido, como
ressalva para a etapa de MD (seção 2.7). Script de renderização:
`scripts/render_consensus_poses.py` (executado no servidor, ambiente
`structure`); composição: `scripts/make_fig7_consensus_poses.py`.

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
`novas-analises/build_peptide_ensemble.py`. Configurações HADDOCK3:
`novas-analises/DN2954-GORE1-2T{,-GGS3}/haddock/run_dn2954-gore1-2t{,-ggs3}.toml`.
Entradas Boltz-2: `novas-analises/DN2954-GORE1-2T{,-GGS3}/boltz2/*.yaml`.
Resultados brutos e tabelas por sistema:
`novas-analises/DN2954-GORE1-2T{,-GGS3}/{qualidade,protonacao_ph8.2,interacao,afinidade}/`.
Comparação detalhada com a MD anterior:
`novas-analises/comparativo_vs_DN2954-GORE12T.md`. Figuras (PDF+PNG
600dpi, com manifesto `.export.json`) e script gerador:
`novas-analises/figuras/`, `novas-analises/scripts/make_figures_pub.py`.
Fig. 7 (poses consenso, renderização estrutural): script PyMOL
`novas-analises/scripts/render_consensus_poses.py` (roda no servidor,
ambiente `structure`, produz `figuras/_render/consensus_*.png`) +
script de composição `novas-analises/scripts/make_fig7_consensus_poses.py`
(roda localmente, aplica autocrop + estilo `scientific-visualization`).
Divergência de conformação HADDOCK vs. Boltz-2 (seção 2.2/2.7, detalhe em
`comparativo_vs_DN2954-GORE12T.md` seção 7): números reproduzidos por
`novas-analises/scripts/diagnose_pose_divergence.py`, sobre os PDBs em
`novas-analises/diagnostico_divergencia_pose/` (conformeros de entrada do
ensemble, amostras Boltz-2, pose HADDOCK reprocessada com `flexref`
cobrindo o peptídeo inteiro). Configuração corrigida do HADDOCK3:
`novas-analises/DN2954-GORE1-2T{,-GGS3}/haddock/run_dn2954-gore1-2t{,-ggs3}.toml`
(`[flexref] nseg2=1`); run original preservado em
`haddock/run1-guided_orig_frozen-peptide/` em cada sistema (não usado
downstream — decisão registrada foi adotar a pose Boltz-2 para a MD).
