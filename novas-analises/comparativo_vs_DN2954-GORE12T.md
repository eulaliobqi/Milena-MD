# Bloco 1 — DN2954 × GORE1-2T / GORE1-2T(GGS)3: validação pré-MD

Numeração de resíduo do receptor sempre na convenção original do PDB
(35–266; His79/Asp126/Ser222/Asp216), a mesma usada em
`assets/samplesheet_dn2954_gore12t.csv`. Poses do Boltz-2 foram remapeadas
(+34) a partir da numeração sequencial 1–232 que ele usa internamente.

## 1. Preparação e reuso de estrutura

- Receptor `DN2954-receptor.pdb` reaproveitado sem refolding; sequência CA
  confere com `DN2954.fa` (faltam só os 3 resíduos C-terminais `TIV`,
  desordenados/não resolvidos).
- `GORE1-2T(GGS)3.fa` é sequência idêntica ao `GORE12T` já usado e
  totalmente analisado em `DN2954-GORE12T/` — usado aqui como o par direto
  de comparação "manteve o padrão?".
- `GORE1-2T` (21 aa, sem linkers) não tinha estrutura prévia — gerado do
  zero (ensemble hélice-α/fita-β/PPII → HADDOCK3; co-folding → Boltz-2).

## 2. Geração de complexo (padrão ouro + fronteira)

| Sistema | Método | Melhor pose | Score/confiança |
|---|---|---|---|
| GORE1-2T | HADDOCK3 (gold) | cluster 1 (10 modelos) | HADDOCK score -67.87 |
| GORE1-2T | Boltz-2 (frontier) | model_0 de 3 | confidence_score 0.90, ipTM 0.81 |
| GORE1-2T-GGS3 | HADDOCK3 (gold) | cluster 1 (10 modelos) | HADDOCK score -56.32 |
| GORE1-2T-GGS3 | Boltz-2 (frontier) | model_0 de 3 | confidence_score 0.87, ipTM 0.73 |

Restrições ambíguas usaram o sítio ativo (His79/Asp126/Ser222/Asp216) como
região ativa do receptor + shell passivo calculado por
`haddock3-restraints passive_from_active` — **esse shell já sobrepõe boa
parte dos resíduos de maior frequência de contato observados na MD anterior
de 100 ns** (ver seção 4), o que valida a escolha de restrição a posteriori.

Confiança do Boltz-2 é consistentemente **mais alta para o peptídeo curto**
(ipTM 0,81/0,83/0,67 nas 3 amostras) do que para o construto ligado GGS3
(ipTM 0,73/0,58/0,46) — esperado, dado que o GGS3 é mais longo e flexível
(75 aa com linkers), mais difícil de "travar" numa única pose confiante.

## 3. Controle de qualidade (nova etapa, antes ausente no Milena-MD)

Geometria da díade catalítica His-NE2···Ser-OG e His-ND1···Asp-OD1/2:
**intacta e consistente em todas as 8 poses avaliadas** (HADDOCK + 3
amostras Boltz-2, para os 2 peptídeos), distâncias 2,7–4,0 Å, dentro da
faixa geometricamente plausível de um par catalítico funcional — nenhuma
pose descartada por tríade rompida.

(Asp216···Ser-OG ficou em ~10–11 Å em todas as poses — isso é **esperado,
não um defeito**: Asp216 é o resíduo de especificidade do fundo do bolsão
S1, que interage com a cadeia lateral do resíduo P1 do substrato/peptídeo,
não faz ponte de H direta com a Ser catalítica.)

## 4. Protonação em pH 8,2 — HIS/ASP/GLU/CYS

- Todo ASP/GLU previsto com pKa 2,0–5,8 → **carregado (padrão)** em pH 8,2,
  em ambos os métodos e ambos os peptídeos — sem exceção.
- Toda CYS (6 por sistema) retornou pKa=99,99 do PROPKA em todas as 4
  poses — sentinela de resíduo **não-titulável / trancado em ponte
  dissulfeto**, consistente com a arquitetura de dissulfetos conhecida da
  tripsina, **sem tiolato livre** em pH 8,2 em nenhuma pose.
- HIS: 5 das 6 His de cada sistema saem neutras (HID/HIE, pKa 4,2–7,0) em
  todas as poses, como esperado.
- **Achado a registrar, não a esconder:** a His79 catalítica sai neutra
  (pKa 7,0–8,0) em 3 das 4 poses — mas na pose HADDOCK do peptídeo curto
  GORE1-2T especificamente, o PROPKA prevê **pKa 8,55, acima do pH de
  referência** (`estado_atribuido = HIP, incomum`). No Boltz-2 (mesmo
  peptídeo) e em ambos os métodos para o GGS3, a mesma His79 sai normal.
  Isso é um efeito local de empacotamento dessa pose específica (não do
  peptídeo em geral) — a reportar como ressalva pontual dessa pose antes de
  usá-la adiante, não como problema do sistema como um todo.

- **Ressalva resolvida — causa raiz identificada por medição direta.** O
  log do PROPKA (`haddock_gore1-2t_propka.log`) atribui o excesso de pKa
  quase inteiramente a um termo de Coulomb com `C- 21 B` (+1,66,
  carboxilato C-terminal do próprio peptídeo, cadeia B), somado ao termo
  já esperado com Asp126 (+1,60/+2,03) — ausente no perfil da mesma His79
  em qualquer outra pose. Medição direta em `Bio.PDB` confirma: His79-NE2
  está a **3,18 Å** do oxigênio C-terminal do resíduo 21 (Arg, última
  posição de GORE1-2T) na pose HADDOCK, contra **13,2–14,1 Å** na pose
  Boltz-2 do mesmo sistema. GORE1-2T é um peptídeo curto (21 aa) sem
  capeamento C-terminal e sem restrição direcional no AIR (tratado como
  totalmente passivo, sem hipótese de qual extremidade encara o sítio) —
  o HADDOCK3 encontrou uma orientação em que a carboxilato livre do C-
  terminal por acaso faz ponte salina com a His catalítica, favorecendo
  sua protonação nessa geometria específica. Isso não se repete no GGS3
  (C-terminal mais distante do bolsão, preso entre os espaçadores) nem em
  nenhuma amostra Boltz-2. **Conclusão: artefato de pose do docking
  rígido+flexível, não uma propriedade real do par receptor-peptídeo** —
  não há razão estrutural para o C-terminal de uma repetição `VLRVLK`
  ocupar esse bolsão em vez do restante da cadeia. **Decisão para a
  próxima etapa:** usar a pose **Boltz-2** (`model_0`, sem esse artefato)
  como estrutura inicial de MD para o sistema GORE1-2T; se a pose HADDOCK
  for usada por algum motivo, His79 deve ser forçada a estado neutro
  (HID/HIE, não HIP) na preparação da topologia, e não o estado literal
  saído do PDB2PQR/PROPKA para essa pose.

## 5. Interação (resíduos de interface, pose estática, corte 4,5 Å)

| Sistema/método | Resíduos de interface (receptor) |
|---|---|
| GORE1-2T · HADDOCK | 56,63,79,84,85,219,222,237-241 |
| GORE1-2T · Boltz-2 | 54,63,64,79,82-84,123,167,170-172,174,197,198,216-219,237-242,247,249 |
| GORE1-2T-GGS3 · HADDOCK | 55,56,63,64,79,80,84,85,120,217-222,238,239,241,242 |
| GORE1-2T-GGS3 · Boltz-2 | 54-56,63,64,79,80,83,84,123,170,171,216-222,236-241,249 |

**Comparação direta com a MD de 100 ns já existente de `DN2954-GORE12T`**
(resíduos de maior frequência de contato: 55,56,57,58,63,99,120-123,
166,167,170,175,197,218,219,220,240,242): o GGS3 novo (mesmo peptídeo)
recupera **8-9 desses resíduos** em cada método (HADDOCK e Boltz-2), e o
GORE1-2T curto (peptídeo novo, sem histórico) recupera **4 resíduos** em
comum com HADDOCK e **9** com Boltz-2 — concentrados na mesma alça
219-222/237-242 ao redor do sítio ativo. **O padrão de interface se
mantém**: as duas ferramentas convergem para a mesma região da MD anterior
já validada, mesmo sem terem visto esse resultado.

## 6. Afinidade (PRODIGY-prot, poses protonadas)

| Sistema | Método | ΔG (kcal/mol) | Kd (M, 25°C) |
|---|---|---|---|
| GORE1-2T | HADDOCK | -6.0 | 3,7×10⁻⁵ |
| GORE1-2T | Boltz-2 | -10.6 | 1,5×10⁻⁸ |
| GORE1-2T-GGS3 | HADDOCK | -6.6 | 1,4×10⁻⁵ |
| GORE1-2T-GGS3 | Boltz-2 | -11.2 | 5,9×10⁻⁹ |

*Ressalva declarada (herdada de `eulalio-pos-doc`): PRODIGY foi calibrado em
interfaces proteína-proteína estruturadas; validade estatística para um
peptídeo curto/não-globular como ligante não está comprovada — reportado
mesmo assim, com a ressalva.*

Duas observações, lado a lado, sem decidir sozinho qual "vence" (mesmo
critério de `consolidate_and_select_pose.py`):
- **Boltz-2 prevê afinidade sistematicamente ~4-5 kcal/mol mais forte que
  HADDOCK**, nos dois peptídeos — provavelmente reflete a maior área de
  interface nas poses de co-folding vs. o refinamento rígido+flexível do
  HADDOCK, não necessariamente maior acurácia.
- **GGS3 (75 aa, com linkers) prevê afinidade levemente mais forte que
  GORE1-2T (21 aa) em ambos os métodos** — direção consistente, mas
  diferença pequena (~0,6 kcal/mol), não decisiva sozinha.

O valor de referência anterior (`DN2954-GORE12T/fe_estimate/free_energy_estimate.txt`,
-104,65 kcal/mol) **não é comparável em valor absoluto** — é uma estimativa
de Interaction Entropy pós-MD, em vácuo, sem termo de solvatação, que o
próprio grupo já marca como não-quantitativo (desvio padrão >10 kT). A
comparação de afinidade válida aqui é de **ranking/direção**, não de
magnitude: nenhum dos dois métodos novos contradiz a existência de uma
interação favorável (ΔG negativo em todos os casos).

## Conclusão do bloco — o padrão se manteve?

**Sim, sem ressalva pendente.** Estrutura, tríade catalítica, interface e
afinidade do peptídeo GGS3 recém-gerado convergem com o que já estava
validado por 100 ns de MD em `DN2954-GORE12T`, usando dois métodos de
docking independentes. O peptídeo novo GORE1-2T (sem linkers) mostra o
mesmo comportamento qualitativo (liga na mesma região do sítio ativo, ΔG
favorável), sem histórico prévio para comparar diretamente. A única
ressalva do bloco — His79 com pKa incomum na pose HADDOCK do GORE1-2T —
foi investigada e resolvida (seção 4): é um artefato de pose (C-terminal
livre do peptídeo curto fazendo ponte salina acidental com a His
catalítica), não uma propriedade real do sistema; decisão tomada de usar
a pose Boltz-2 para esse sistema em MD (ou forçar His79 neutra se a pose
HADDOCK for usada).

**Gate para o próximo bloco (MD triplicata):** revisão humana deste
documento pelo usuário.
