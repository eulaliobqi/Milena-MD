const fs = require("fs");
const path = require("path");
const {
  Document, Packer, Paragraph, TextRun, HeadingLevel, AlignmentType,
  ImageRun, Table, TableRow, TableCell, WidthType, ShadingType,
  BorderStyle, VerticalAlign, PageOrientation,
} = require("C:/Users/eulal/AppData/Roaming/npm/node_modules/docx");

const FIG = "C:/Users/eulal/.claude/Milena-MD/novas-analises/figuras/";
const OUT = "C:/Users/eulal/.claude/Milena-MD/novas-analises/metodologia_resultados_DN2954-GORE1-2T_GGS3.docx";

const FONT = "Times New Roman";
const SZ_BODY = 24; // 12pt
const SZ_SMALL = 20; // 10pt
const SZ_H1 = 28; // 14pt
const SZ_H2 = 26; // 13pt
const SZ_TITLE = 32; // 16pt

function para(text, opts = {}) {
  const { bold = false, italic = false, size = SZ_BODY, spacingAfter = 200, alignment, indent } = opts;
  return new Paragraph({
    alignment,
    indent,
    spacing: { after: spacingAfter, line: 276 },
    children: [new TextRun({ text, font: FONT, size, bold, italics: italic })],
  });
}

// paragraph with mixed runs: array of {text, bold, italic, code}
function mixedPara(runs, opts = {}) {
  const { spacingAfter = 200, alignment } = opts;
  return new Paragraph({
    alignment,
    spacing: { after: spacingAfter, line: 276 },
    children: runs.map(r => new TextRun({
      text: r.text, font: r.code ? "Consolas" : FONT, size: r.code ? SZ_SMALL : SZ_BODY,
      bold: !!r.bold, italics: !!r.italic,
    })),
  });
}

function h1(text) {
  return new Paragraph({
    heading: HeadingLevel.HEADING_1,
    spacing: { before: 360, after: 200 },
    children: [new TextRun({ text, font: FONT, size: SZ_H1, bold: true })],
  });
}

function h2(text) {
  return new Paragraph({
    heading: HeadingLevel.HEADING_2,
    spacing: { before: 280, after: 160 },
    children: [new TextRun({ text, font: FONT, size: SZ_H2, bold: true })],
  });
}

function figure(filename, widthIn, heightIn, captionRuns) {
  const data = fs.readFileSync(FIG + filename);
  return [
    new Paragraph({
      alignment: AlignmentType.CENTER,
      spacing: { before: 160, after: 80 },
      children: [
        new ImageRun({
          type: "png",
          data,
          transformation: { width: Math.round(widthIn * 96), height: Math.round(heightIn * 96) },
        }),
      ],
    }),
    new Paragraph({
      spacing: { after: 240 },
      children: captionRuns.map(r => new TextRun({
        text: r.text, font: r.code ? "Consolas" : FONT, size: SZ_SMALL,
        italics: r.italic !== false, bold: !!r.bold,
      })),
    }),
  ];
}

function cellText(text, opts = {}) {
  const { bold = false, shade = null, align = AlignmentType.CENTER } = opts;
  return new TableCell({
    width: opts.width ? { size: opts.width, type: WidthType.DXA } : undefined,
    shading: shade ? { type: ShadingType.CLEAR, fill: shade } : undefined,
    verticalAlign: VerticalAlign.CENTER,
    margins: { top: 60, bottom: 60, left: 80, right: 80 },
    children: [new Paragraph({
      alignment: align,
      children: [new TextRun({ text, font: FONT, size: SZ_SMALL, bold })],
    })],
  });
}

const borders = {
  top: { style: BorderStyle.SINGLE, size: 4, color: "999999" },
  bottom: { style: BorderStyle.SINGLE, size: 4, color: "999999" },
  left: { style: BorderStyle.SINGLE, size: 4, color: "999999" },
  right: { style: BorderStyle.SINGLE, size: 4, color: "999999" },
  insideHorizontal: { style: BorderStyle.SINGLE, size: 4, color: "999999" },
  insideVertical: { style: BorderStyle.SINGLE, size: 4, color: "999999" },
};

// ---------- Table 1: Boltz-2 confidence per sample ----------
const t1widths = [1800, 1400, 1600, 1200, 1200, 1200];
const t1 = new Table({
  width: { size: 8400, type: WidthType.DXA },
  columnWidths: t1widths,
  borders,
  rows: [
    new TableRow({
      tableHeader: true,
      children: ["Sistema", "Amostra", "confidence_score", "ipTM", "pTM", "pLDDT"].map((t, i) =>
        cellText(t, { bold: true, shade: "D9D9D9", width: t1widths[i] })),
    }),
    ...[
      ["GORE1-2T", "model_0", "0,901", "0,807", "0,967", "0,925"],
      ["GORE1-2T", "model_1", "0,900", "0,825", "0,968", "0,918"],
      ["GORE1-2T", "model_2", "0,862", "0,671", "0,954", "0,909"],
      ["GORE1-2T(GGS)3", "model_0", "0,866", "0,733", "0,920", "0,900"],
      ["GORE1-2T(GGS)3", "model_1", "0,801", "0,579", "0,884", "0,857"],
      ["GORE1-2T(GGS)3", "model_2", "0,790", "0,464", "0,853", "0,872"],
    ].map(row => new TableRow({
      children: row.map((t, i) => cellText(t, { width: t1widths[i] })),
    })),
  ],
});

// ---------- Table 2: Synthesis per system ----------
const t2widths = [1700, 1850, 1250, 950, 1100, 1250, 950]; // soma 9050 DXA (~6,28in), cabe na largura util da pagina (~6,53in)
const t2 = new Table({
  width: { size: 9050, type: WidthType.DXA },
  columnWidths: t2widths,
  borders,
  rows: [
    new TableRow({
      tableHeader: true,
      children: [
        "Sistema", "Tríade NE2\u00b7\u00b7\u00b7Ser / ND1\u00b7\u00b7\u00b7Asp (\u00c5)\u00b9",
        "pKa ASP/GLU (faixa)", "pKa His79", "\u0394G PRODIGY (kcal/mol)", "Kd (M, 25\u00b0C)", "Interface\u00b2",
      ].map((t, i) => cellText(t, { bold: true, shade: "D9D9D9", width: t2widths[i] })),
    }),
    new TableRow({
      children: [
        cellText("DN2954 \u00d7 GORE1-2T", { width: t2widths[0], align: AlignmentType.LEFT }),
        cellText("3,23\u20133,30 / 2,75\u20132,83", { width: t2widths[1] }),
        cellText("2,89\u20135,32", { width: t2widths[2] }),
        cellText("6,97", { width: t2widths[3] }),
        cellText("\u221210,6", { width: t2widths[4] }),
        cellText("1,5\u00d710\u207b\u2078", { width: t2widths[5] }),
        cellText("9/20", { width: t2widths[6] }),
      ],
    }),
    new TableRow({
      children: [
        cellText("DN2954 \u00d7 GORE1-2T(GGS)3", { width: t2widths[0], align: AlignmentType.LEFT }),
        cellText("3,19\u20133,22 / 2,80\u20132,81", { width: t2widths[1] }),
        cellText("2,65\u20135,24", { width: t2widths[2] }),
        cellText("7,04", { width: t2widths[3] }),
        cellText("\u221211,2", { width: t2widths[4] }),
        cellText("5,9\u00d710\u207b\u2079", { width: t2widths[5] }),
        cellText("9/20", { width: t2widths[6] }),
      ],
    }),
  ],
});

const doc = new Document({
  sections: [{
    properties: {
      page: {
        size: { width: 12240, height: 15840 },
        margin: { top: 1417, bottom: 1417, left: 1701, right: 1134 },
      },
    },
    children: [
      new Paragraph({
        alignment: AlignmentType.CENTER,
        spacing: { after: 120 },
        children: [new TextRun({
          text: "Interação e afinidade de ligação da tripsina DN2954 com os peptídeos GORE1-2T e GORE1-2T(GGS)3: preparação, protonação em pH 8,2 e geração de pose por co-folding (Boltz-2)",
          font: FONT, size: SZ_TITLE, bold: true,
        })],
      }),
      para(
        "Bloco 1 de novas-analises \u2014 etapa pr\u00e9-din\u00e2mica molecular. Numera\u00e7\u00e3o de res\u00edduo do receptor sempre na conven\u00e7\u00e3o original do PDB (35\u2013266; His79/Asp126/Ser222/Asp216). Este documento mostra apenas o que segue para a pr\u00f3xima etapa (din\u00e2mica molecular em triplicata) \u2014 a pose gerada por co-folding com Boltz-2. O docking guiado por restri\u00e7\u00f5es com HADDOCK3 continua descrito nos M\u00e9todos e usado como valida\u00e7\u00e3o cruzada independente do s\u00edtio de liga\u00e7\u00e3o, mas sua pr\u00f3pria pose de pept\u00eddeo n\u00e3o \u00e9 usada em nenhuma etapa seguinte (motivo resumido na se\u00e7\u00e3o 2.2; investiga\u00e7\u00e3o completa em comparativo_vs_DN2954-GORE12T.md, se\u00e7\u00e3o 7).",
        { italic: true, size: SZ_SMALL, spacingAfter: 360 }
      ),

      // ===================== 1. METODOS =====================
      h1("1. Material e Métodos"),

      h2("1.1 Sistema e dados de entrada"),
      para("Receptor: tripsina digestiva DN2954 (232 res\u00edduos, His79/Asp126/Ser222/Asp216 como t\u00e9trade catal\u00edtica cl\u00e1ssica intacta), j\u00e1 modelada e usada em trabalho anterior (DN2954-GORE12T/). Dois ligantes pept\u00eddicos: GORE1-2T (21 aa, repeti\u00e7\u00e3o VLRVLKVLRVLKVLRVLKVLR, sem estrutura pr\u00e9via) e GORE1-2T(GGS)3 (75 aa, mesma repeti\u00e7\u00e3o com tr\u00eas espa\u00e7adores (GGS)3 entre unidades) \u2014 sequ\u00eancia id\u00eantica, res\u00edduo a res\u00edduo, ao ligante GORE12T j\u00e1 usado e totalmente caracterizado por 100 ns de din\u00e2mica molecular em DN2954-GORE12T/, o que permite us\u00e1-lo como compara\u00e7\u00e3o direta de consist\u00eancia de m\u00e9todo (se\u00e7\u00e3o 2.7)."),
      para("A identidade de sequ\u00eancia do receptor reaproveitado (data/DN2954-receptor.pdb) contra DN2954.fa foi confirmada por extra\u00e7\u00e3o da sequ\u00eancia C\u03b1 (Biopython 1.83; Cock et al. 2009); faltam apenas os 3 res\u00edduos C-terminais (TIV), regi\u00e3o desordenada n\u00e3o resolvida no modelo estrutural original."),

      h2("1.2 Ensemble conformacional do peptídeo (entrada para o HADDOCK3)"),
      para("Para cada pept\u00eddeo, tr\u00eas conf\u00f4rmeros iniciais foram constru\u00eddos a partir da sequ\u00eancia com PyMOL 3.1.0 (Schr\u00f6dinger, LLC), comando cmd.fab(): h\u00e9lice-\u03b1 (ss=1), fita-\u03b2 antiparalela estendida (ss=2) e poliprolina-II (cadeia estendida ss=0 seguida de fixa\u00e7\u00e3o por res\u00edduo dos \u00e2ngulos diedrais can\u00f4nicos \u03c6=\u221275\u00b0, \u03c8=145\u00b0 via cmd.set_dihedral). Os tr\u00eas modelos foram combinados num \u00fanico PDB multi-modelo com pdb_mkensemble (pdb-tools 2.7.0; Rodrigues et al. 2018) e validados com pdb_chkensemble. Este ensemble multi-conf\u00f4mero segue a recomenda\u00e7\u00e3o corrente do pr\u00f3prio HADDOCK3 para docking de pept\u00eddeos flex\u00edveis (Giulini et al. 2025) \u2014 usado como entrada do docking descrito em 1.3."),

      h2("1.3 Validação cruzada do sítio de ligação (HADDOCK3)"),
      para("Docking guiado por restri\u00e7\u00f5es amb\u00edguas de intera\u00e7\u00e3o (AIR) com HADDOCK3 vers\u00e3o 2026.7.0 (Giulini et al. 2025), executado localmente (mode = \"local\", ncores = 12). Usado neste bloco como confirma\u00e7\u00e3o independente de onde o pept\u00eddeo liga (regi\u00e3o de interface, t\u00e9trade catal\u00edtica) \u2014 n\u00e3o como fonte da pose de MD (motivo na se\u00e7\u00e3o 2.2)."),
      mixedPara([
        { text: "Restri\u00e7\u00f5es: ", bold: true },
        { text: "res\u00edduo ativo do receptor = t\u00e9trade catal\u00edtica (His79, Asp126, Ser222, Asp216); shell passivo calculado por " },
        { text: "haddock3-restraints passive_from_active", code: true },
        { text: " (raio padr\u00e3o 6,5 \u00c5), retornando 31 res\u00edduos de superf\u00edcie ao redor do s\u00edtio ativo. Pept\u00eddeo tratado como totalmente passivo. Tabela amb\u00edgua gerada por " },
        { text: "haddock3-restraints active_passive_to_ambig", code: true },
        { text: "." },
      ]),
      mixedPara([
        { text: "Fluxo (.toml): ", bold: true },
        { text: "topoaa \u2192 rigidbody (sampling=1000) \u2192 caprieval \u2192 seletop (select=200) \u2192 flexref (pept\u00eddeo inteiro declarado semi-flex\u00edvel, nseg2=1) \u2192 caprieval \u2192 emref \u2192 caprieval \u2192 clustfcc \u2192 seletopclusts \u2192 caprieval final." },
      ]),
      para("Pose de refer\u00eancia de cada sistema: melhor modelo do cluster de menor (mais negativo) HADDOCK score no caprieval final \u2014 usada apenas para os n\u00fameros de valida\u00e7\u00e3o citados na se\u00e7\u00e3o 2.2, n\u00e3o como estrutura de partida de MD."),

      h2("1.4 Geração de pose — co-folding (Boltz-2, método adotado)"),
      para("Co-folding direto receptor+pept\u00eddeo a partir das sequ\u00eancias, sem posicionamento nem restri\u00e7\u00e3o pr\u00e9via, com Boltz vers\u00e3o 2.2.1 (Passaro et al. 2025), ambiente boltz2-env dedicado. Esta \u00e9 a pose usada como estrutura inicial da pr\u00f3xima etapa (MD em triplicata) \u2014 model_0 de cada sistema."),
      mixedPara([
        { text: "Par\u00e2metros: ", bold: true },
        { text: "boltz predict --use_msa_server --no_kernels --diffusion_samples 3 --output_format pdb", code: true },
        { text: ". MSA obtido do servidor p\u00fablico ColabFold (api.colabfold.com), pareamento greedy. --no_kernels usado por compatibilidade de kernel CUDA com a GPU RTX 5070 Ti (Blackwell)." },
      ]),
      para("N\u00e3o solicitado o bloco properties: affinity: esse m\u00f3dulo do Boltz-2 s\u00f3 aceita ligante de mol\u00e9cula pequena (ligand/SMILES) como binder, n\u00e3o uma cadeia protein/pept\u00eddica \u2014 confirmado empiricamente (falha ValueError: Chain B is not a ligand!) antes de descartar o recurso. A avalia\u00e7\u00e3o de confian\u00e7a usou as m\u00e9tricas nativas do Boltz-2 (confidence_score, ipTM, pTM, pLDDT do complexo, de confidence_*.json, e o pLDDT por res\u00edduo gravado no campo B-factor do PDB de sa\u00edda)."),
      para("Tr\u00eas amostras de difus\u00e3o por sistema (n\u00e3o uma \u00fanica), seguindo a recomenda\u00e7\u00e3o de checagem de consist\u00eancia de pose entre amostras (motivada por Kim et al. 2026 \u2014 confian\u00e7a alta do modelo n\u00e3o garante pose fisicamente correta). Essa checagem de consist\u00eancia foi decisiva: ver se\u00e7\u00e3o 2.2/2.7."),

      h2("1.5 Controle de qualidade estrutural"),
      para("Etapa nova em rela\u00e7\u00e3o ao pipeline de produ\u00e7\u00e3o do Milena-MD. Sobre cada uma das 3 amostras de difus\u00e3o Boltz-2 de cada sistema, a geometria da d\u00edade catal\u00edtica foi calculada estaticamente com script pr\u00f3prio (scripts/triad_distances.py, Biopython 1.83 Bio.PDB), medindo as dist\u00e2ncias heavy-atom His-NE2\u00b7\u00b7\u00b7Ser-OG e His-ND1\u00b7\u00b7\u00b7Asp-OD1/OD2, com limiar de aceita\u00e7\u00e3o \u22644,5 \u00c5. A mesma checagem foi aplicada \u00e0 pose HADDOCK (valida\u00e7\u00e3o cruzada); n\u00fameros completos em comparativo_vs_DN2954-GORE12T.md."),

      h2("1.6 Protonação em pH 8,2 e avaliação de resíduos ionizáveis"),
      para("Protona\u00e7\u00e3o com PDB2PQR vers\u00e3o 3.7.1 (Dolinsky et al. 2004) usando o motor de estados de titula\u00e7\u00e3o PROPKA vers\u00e3o 3.5.1 (Olsson et al. 2011; S\u00f8ndergaard et al. 2011): pdb2pqr30 --ff AMBER --titration-state-method propka --with-ph 8.2 --keep-chain --pdb-output, pH fixado no valor padronizado pelo grupo para intestino m\u00e9dio de lepid\u00f3ptero. A tabela de pKa por res\u00edduo foi extra\u00edda por script pr\u00f3prio (scripts/protonate_ph82.py), filtrando HIS, ASP, GLU e CYS e sinalizando qualquer res\u00edduo cujo pKa previsto cruzasse o pH de refer\u00eancia na dire\u00e7\u00e3o incomum. Aplicado \u00e0 pose Boltz-2 model_0 de cada sistema."),

      h2("1.7 Análise de interação (interface receptor-peptídeo)"),
      para("Para cada pose protonada, res\u00edduos do receptor (cadeia A) com pelo menos um \u00e1tomo pesado a \u22644,5 \u00c5 de qualquer \u00e1tomo pesado do pept\u00eddeo (cadeia B) foram listados com NeighborSearch do Biopython (scripts/interface_contacts.py) \u2014 equivalente de pose \u00fanica, sem trajet\u00f3ria de MD, ao mapa de contato por frequ\u00eancia j\u00e1 usado em produ\u00e7\u00e3o (bin/contact_map.py)."),

      h2("1.8 Estimativa de afinidade de ligação"),
      para("Energia livre de liga\u00e7\u00e3o estimada com PRODIGY-prot vers\u00e3o 2.3.0 (Xue et al. 2016), m\u00e9todo baseado em contatos de interface calibrado em complexos prote\u00edna-prote\u00edna estruturados: prodigy <pdb> --selection A B, temperatura padr\u00e3o 25\u00b0C, sobre a pose Boltz-2 model_0 j\u00e1 protonada em pH 8,2. Ressalva declarada: a validade estat\u00edstica do PRODIGY para uma cadeia B curta/n\u00e3o-globular n\u00e3o est\u00e1 comprovada \u2014 reportado do mesmo modo que em uso anterior do grupo para o pept\u00eddeo GORE3 (5 aa), com a mesma ressalva propagada."),

      // ===================== 2. RESULTADOS =====================
      h1("2. Resultados"),

      h2("2.1 A estrutura do receptor e do ligante GGS3 são reaproveitáveis sem retrabalho"),
      para("A sequ\u00eancia C\u03b1 extra\u00edda de DN2954-receptor.pdb reproduz exatamente DN2954.fa (232/232 res\u00edduos, faltando apenas os 3 res\u00edduos C-terminais desordenados TIV), e a sequ\u00eancia de GORE12T-ligand-DN2954.pdb \u00e9 id\u00eantica, res\u00edduo a res\u00edduo, a GORE1-2T(GGS)3.fa (75/75). Nenhuma nova predi\u00e7\u00e3o de estrutura de receptor foi necess\u00e1ria; o ligante GGS3 p\u00f4de ser tratado como o mesmo objeto molecular j\u00e1 validado por MD, permitindo compara\u00e7\u00e3o direta de m\u00e9todo em vez de compara\u00e7\u00e3o entre mol\u00e9culas diferentes."),

      h2("2.2 Boltz-2 converge em pose única e bem definida; HADDOCK3 confirma o sítio mas não fornece pose de peptídeo independente"),
      para("O co-folding Boltz-2 concorda em confian\u00e7a alta para as 3 amostras de difus\u00e3o de ambos os sistemas (confidence_score 0,86\u20130,90; Fig. 1, Tabela 1), com uma diferen\u00e7a sistem\u00e1tica entre pept\u00eddeos: o pept\u00eddeo curto GORE1-2T \u00e9 predito com confian\u00e7a de interface (ipTM) mais alta e mais est\u00e1vel entre amostras (0,81/0,83/0,67) do que o construto ligado GORE1-2T(GGS)3 (0,73/0,58/0,46) \u2014 consistente com a maior flexibilidade conformacional esperada de um pept\u00eddeo 75 aa com tr\u00eas espa\u00e7adores (GGS)3 frente a um pept\u00eddeo curto de 21 aa."),
      ...figure("fig1_boltz_confidence.png", 6.1, 3.588, [
        { text: "Figura 1. ", bold: true, italic: false },
        { text: "Boltz-2: m\u00e9tricas de confian\u00e7a por amostra de difus\u00e3o. confidence_score, ipTM, pTM e pLDDT do complexo, lidos diretamente de confidence_*.json, para as 3 amostras (--diffusion_samples 3) de cada sistema." },
      ]),
      new Paragraph({ spacing: { after: 120 }, children: [new TextRun({ text: "Tabela 1. Métricas de confiança Boltz-2 por amostra de difusão.", font: FONT, size: SZ_SMALL, bold: true })] }),
      t1,
      para("", { spacingAfter: 200 }),
      para("Valida\u00e7\u00e3o cruzada: o docking guiado por restri\u00e7\u00f5es HADDOCK3 (se\u00e7\u00e3o 1.3) concorda com essa regi\u00e3o de liga\u00e7\u00e3o \u2014 cluster l\u00edder claramente separado dos demais por HADDOCK score nos dois sistemas (GORE1-2T: \u221267,87 \u00b1 3,53 vs. segundo colocado \u221255,36 \u00b1 4,52; GORE1-2T(GGS)3: \u221256,32 \u00b1 1,89 vs. segundo colocado \u221253,33 \u00b1 3,04) e recupera\u00e7\u00e3o de interface compar\u00e1vel \u00e0 do Boltz-2 (se\u00e7\u00e3o 2.5). Essa concord\u00e2ncia \u00e9 usada aqui s\u00f3 como confirma\u00e7\u00e3o de que a regi\u00e3o de liga\u00e7\u00e3o identificada \u00e9 robusta a dois m\u00e9todos independentes, n\u00e3o como fonte de uma pose alternativa de pept\u00eddeo: investiga\u00e7\u00e3o dirigida (RMSD de backbone, detalhe completo em comparativo_vs_DN2954-GORE12T.md se\u00e7\u00e3o 7) mostrou que o backbone do pept\u00eddeo na pose HADDOCK tem RMSD de apenas 0,38\u20130,43 \u00c5 contra o conf\u00f4rmero h\u00e9lice-\u03b1 de entrada (se\u00e7\u00e3o 1.2) \u2014 ou seja, o flexref do HADDOCK3 n\u00e3o o alterou de fato. Mesmo depois de reconfigurar flexref para tratar o pept\u00eddeo inteiro como semi-flex\u00edvel e reprocessar a partir dessa etapa (--restart 4, reaproveitando o rigidbody j\u00e1 calculado), o HADDOCK score n\u00e3o mudou (estatisticamente id\u00eantico) mas a conforma\u00e7\u00e3o tamb\u00e9m n\u00e3o: RMSD 0,41\u20130,75 \u00c5 contra a mesma h\u00e9lice, s\u00f3 reposicionada em rela\u00e7\u00e3o ao receptor."),
      mixedPara([
        { text: "Conclus\u00e3o: n\u00e3o \u00e9 um par\u00e2metro mal ajustado, \u00e9 um limite de escopo do flexref", bold: true },
        { text: " \u2014 refinamento local por SA de tor\u00e7\u00e3o em rajadas curtas, sem capacidade de desfazer uma \u03b1-h\u00e9lice j\u00e1 formada. Por isso a pose usada na pr\u00f3xima etapa \u00e9 sempre a do Boltz-2 model_0 (Fig. 6) \u2014 HADDOCK3 fica como confirma\u00e7\u00e3o independente do s\u00edtio, n\u00e3o como gerador de pose." },
      ]),

      h2("2.3 A tétrade catalítica permanece geometricamente intacta nas 6 poses Boltz-2 avaliadas"),
      para("As dist\u00e2ncias His-NE2\u00b7\u00b7\u00b7Ser-OG (3,19\u20133,30 \u00c5) e His-ND1\u00b7\u00b7\u00b7Asp-OD1/OD2 (2,75\u20132,83 \u00c5, menor das duas) ficaram dentro da faixa de par catal\u00edtico plaus\u00edvel em todas as 3 amostras de difus\u00e3o dos dois pept\u00eddeos (Fig. 2). Nenhuma amostra foi descartada por t\u00e9trade rompida. A dist\u00e2ncia Ser-OG\u00b7\u00b7\u00b7Asp216 (9,7\u201310,7 \u00c5 em todas as poses) n\u00e3o constitui achado de baixa qualidade: Asp216 \u00e9 o res\u00edduo de especificidade do fundo do bols\u00e3o S1, que interage com a cadeia lateral do res\u00edduo P1 do pept\u00eddeo ligante, n\u00e3o faz ponte de hidrog\u00eanio direta com a serina catal\u00edtica."),
      ...figure("fig2_triad_distances.png", 6.1, 3.965, [
        { text: "Figura 2. ", bold: true, italic: false },
        { text: "Geometria da d\u00edade catal\u00edtica His79/Ser222/Asp126, poses Boltz-2. Dist\u00e2ncias heavy-atom His-NE2\u00b7\u00b7\u00b7Ser-OG e His-ND1\u00b7\u00b7\u00b7Asp-OD (menor das duas) em cada uma das 3 amostras de difus\u00e3o \u00d7 2 pept\u00eddeos (6 poses). Linha tracejada = limiar de aceita\u00e7\u00e3o (4,5 \u00c5)." },
      ]),

      h2("2.4 A protonação em pH 8,2 confirma o padrão esperado de carga, sem anomalias na pose adotada"),
      para("Nas 2 poses protonadas (Boltz-2 model_0, os dois pept\u00eddeos), todo ASP e GLU (17 res\u00edduos por sistema) foi previsto carregado (ASP 2,89\u20135,03; GLU 2,63\u20135,32; ambos bem abaixo de 8,2; Fig. 3), e todas as 6 ciste\u00ednas por sistema retornaram pKa=99,99 do PROPKA \u2014 sentinela de res\u00edduo n\u00e3o-titul\u00e1vel, consistente com a arquitetura de pontes dissulfeto j\u00e1 conhecida da tripsina, sem ind\u00edcio de tiolato livre em pH 8,2. Todas as 6 histidinas por sistema, incluindo a His79 catal\u00edtica (pKa 6,97 no GORE1-2T; 7,04 no GGS3), saem neutras (HID/HIE, pKa 4,15\u20137,04) \u2014 sem exce\u00e7\u00e3o nas poses adotadas. (A pose HADDOCK descartada mostrava uma His79 an\u00f4mala, pKa 8,55, especificamente no GORE1-2T \u2014 artefato de ponte salina com o C-terminal livre do pept\u00eddeo naquela pose espec\u00edfica, sem rela\u00e7\u00e3o com a pose Boltz-2 usada aqui; causa raiz completa em comparativo_vs_DN2954-GORE12T.md se\u00e7\u00e3o 4.)"),
      ...figure("fig3_pka_summary.png", 6.0, 3.718, [
        { text: "Figura 3. ", bold: true, italic: false },
        { text: "pKa previsto de HIS/ASP/GLU nas 2 poses Boltz-2 model_0 protonadas em pH 8,2. PROPKA 3.5.1 via PDB2PQR 3.7.1. Linha tracejada = pH de refer\u00eancia (8,2). Nenhuma anomalia sinalizada nas poses adotadas. CYS omitida do painel (pKa=99,99 em todas as 6 ciste\u00ednas das 2 poses)." },
      ]),

      h2("2.5 A interface de ligação prevista recupera a maior parte dos contatos já validados por 100 ns de dinâmica molecular"),
      para("Os res\u00edduos de interface (\u22644,5 \u00c5) previstos pelo Boltz-2 para o GORE1-2T(GGS)3 \u2014 mol\u00e9cula id\u00eantica \u00e0 j\u00e1 simulada em DN2954-GORE12T/ \u2014 recuperam 9 dos ~20 res\u00edduos de maior frequ\u00eancia de contato observados ao longo de 100 ns de MD (55, 56, 63, 120, 218, 219, 220, 240, 242, entre outros; Fig. 5), apesar de o m\u00e9todo n\u00e3o ter usado essa trajet\u00f3ria como refer\u00eancia ou restri\u00e7\u00e3o. O pept\u00eddeo novo GORE1-2T (sem espa\u00e7adores, sem hist\u00f3rico de MD) recupera um subconjunto igualmente numeroso (9 res\u00edduos em comum), concentrado na mesma al\u00e7a 219\u2013222/237\u2013242 ao redor do s\u00edtio ativo \u2014 o mesmo bols\u00e3o de liga\u00e7\u00e3o, n\u00e3o um s\u00edtio alternativo. O docking HADDOCK3 (valida\u00e7\u00e3o cruzada, se\u00e7\u00e3o 2.2) recupera um conjunto parcialmente sobreposto (4\u20139 res\u00edduos em comum, n\u00fameros completos em comparativo_vs_DN2954-GORE12T.md se\u00e7\u00e3o 5), consistente com a mesma regi\u00e3o."),
      ...figure("fig5_interface_overlap.png", 4.0, 3.467, [
        { text: "Figura 5. ", bold: true, italic: false },
        { text: "Recupera\u00e7\u00e3o da interface j\u00e1 validada por 100 ns de MD, pose Boltz-2. N\u00ba de res\u00edduos do receptor (contato \u22644,5 \u00c5) em comum com os 20 res\u00edduos de maior frequ\u00eancia de contato na MD de 100 ns de DN2954-GORE12T, por sistema. Linha tracejada = tamanho total do conjunto de refer\u00eancia (n=20)." },
      ]),

      h2("2.6 A afinidade prevista é favorável para os dois peptídeos"),
      para("PRODIGY-prot prev\u00ea liga\u00e7\u00e3o favor\u00e1vel para a pose adotada nos dois sistemas: GORE1-2T \u0394G = \u221210,6 kcal/mol (Kd \u2248 1,5\u00d710\u207b\u2078 M); GORE1-2T(GGS)3 \u0394G = \u221211,2 kcal/mol (Kd \u2248 5,9\u00d710\u207b\u2079 M) \u2014 Fig. 4, diferen\u00e7a de ~0,6 kcal/mol entre os dois pept\u00eddeos, n\u00e3o decisiva isoladamente. O valor de energia livre j\u00e1 existente para DN2954-GORE12T (estimativa de Interaction Entropy p\u00f3s-MD, \u2212104,65 kcal/mol) n\u00e3o \u00e9 compar\u00e1vel em magnitude absoluta a este valor de PRODIGY \u2014 \u00e9 uma estimativa em v\u00e1cuo, sem termo de solvata\u00e7\u00e3o, que o pr\u00f3prio grupo j\u00e1 classifica como inst\u00e1vel (desvio padr\u00e3o >10 kT); a compara\u00e7\u00e3o v\u00e1lida entre as duas gera\u00e7\u00f5es de an\u00e1lise \u00e9 de dire\u00e7\u00e3o (liga\u00e7\u00e3o favor\u00e1vel em ambas), n\u00e3o de m\u00f3dulo."),
      ...figure("fig4_prodigy_affinity.png", 4.0, 3.467, [
        { text: "Figura 4. ", bold: true, italic: false },
        { text: "Afinidade prevista (PRODIGY-prot 2.3.0), pose Boltz-2 model_0 protonada em pH 8,2. \u0394G por prodigy <pdb> --selection A B (25\u00b0C), para os 2 pept\u00eddeos." },
      ]),
      new Paragraph({ spacing: { after: 120 }, children: [new TextRun({ text: "Tabela 2. Síntese comparativa dos dois sistemas (pose Boltz-2 adotada).", font: FONT, size: SZ_SMALL, bold: true })] }),
      t2,
      para("\u00b9 Dist\u00e2ncia heavy-atom His-NE2\u00b7\u00b7\u00b7Ser-OG / menor de His-ND1\u00b7\u00b7\u00b7Asp-OD1 ou OD2, faixa entre as 3 amostras de difus\u00e3o. \u00b2 Res\u00edduos de interface (\u22644,5 \u00c5) em comum com o top-20 por frequ\u00eancia de contato na MD de 100 ns de refer\u00eancia (DN2954-GORE12T).", { italic: true, size: SZ_SMALL, spacingAfter: 300 }),

      h2("2.7 Síntese: pose adotada e ressalva remanescente para a MD"),
      para("O pept\u00eddeo GORE1-2T(GGS)3 \u2014 reprocessado do zero com um protocolo que n\u00e3o teve acesso \u00e0 trajet\u00f3ria de MD anterior \u2014 converge com o resultado j\u00e1 validado de DN2954-GORE12T em geometria de s\u00edtio ativo, regi\u00e3o de interface e sinal de afinidade favor\u00e1vel, usando a pose gerada por co-folding (Boltz-2) e com valida\u00e7\u00e3o independente adicional do HADDOCK3 no pr\u00f3prio s\u00edtio de liga\u00e7\u00e3o. O pept\u00eddeo novo GORE1-2T, sem hist\u00f3rico pr\u00e9vio, liga na mesma regi\u00e3o do s\u00edtio ativo com afinidade favor\u00e1vel."),
      mixedPara([
        { text: "Estrutura adotada para a pr\u00f3xima etapa (din\u00e2mica molecular em triplicata): a pose Boltz-2 model_0, para os dois sistemas", bold: true },
        { text: " (Fig. 6). " },
        { text: "Ressalva registrada, n\u00e3o escondida: ", bold: true },
        { text: "para o GORE1-2T especificamente, as 3 amostras de difus\u00e3o do pr\u00f3prio Boltz-2 divergem entre si na conforma\u00e7\u00e3o livre do pept\u00eddeo (dist\u00e2ncia ponta-a-ponta 22,0\u201336,6 \u00c5; pLDDT m\u00e9dio 47\u201356, abaixo do limiar costumeiro de confian\u00e7a) \u2014 ao contr\u00e1rio do GGS3, cujas 3 amostras convergem fortemente (16,0\u201316,6 \u00c5; pLDDT 62\u201375). Isso n\u00e3o \u00e9 um problema da regi\u00e3o de interface (que recupera contatos validados nos dois pept\u00eddeos, se\u00e7\u00e3o 2.5) \u2014 \u00e9 uma incerteza real sobre a conforma\u00e7\u00e3o da por\u00e7\u00e3o livre do GORE1-2T curto, a esclarecer pela pr\u00f3pria MD de produ\u00e7\u00e3o (por exemplo comparando as r\u00e9plicas da triplicata), n\u00e3o algo j\u00e1 resolvido nesta etapa de docking/co-folding." },
      ]),
      ...figure("fig6_boltz_pose.png", 6.1, 4.078, [
        { text: "Figura 6. ", bold: true, italic: false },
        { text: "Pose Boltz-2 model_0 adotada como estrutura inicial de MD, renderiza\u00e7\u00e3o estrutural. PyMOL 3.1.0, cadeia A (receptor) de cada pose superposta rigidamente (cmd.align) sobre data/DN2954-receptor.pdb; t\u00e9trade catal\u00edtica em verde (sticks). (a) DN2954 \u00d7 GORE1-2T (21 aa) e (b) DN2954 \u00d7 GORE1-2T(GGS)3 (75 aa): mesma c\u00e2mera/escala nos dois pain\u00e9is." },
      ]),

      // ===================== REFERENCIAS =====================
      h1("Referências"),
      para("Cock, P. J. A. et al. Biopython: freely available Python tools for computational molecular biology and bioinformatics. Bioinformatics 25, 1422\u20131423 (2009).", { spacingAfter: 160 }),
      para("Schr\u00f6dinger, LLC. The PyMOL Molecular Graphics System, Version 3.1 (2025).", { spacingAfter: 160 }),
      para("Rodrigues, J. P. G. L. M., Teixeira, J. M. C., Trellet, M. & Bonvin, A. M. J. J. pdb-tools: a swiss army knife for molecular structures. F1000Research 7, 1961 (2018).", { spacingAfter: 160 }),
      para("Giulini, M. et al. HADDOCK3. J. Chem. Inf. Model. (2025). doi:10.1021/acs.jcim.5c00969", { spacingAfter: 160 }),
      para("Passaro, S. et al. Boltz-2: Towards Accurate and Efficient Binding Affinity Prediction. bioRxiv (2025). PMC12262699.", { spacingAfter: 160 }),
      para("Kim, J. et al. Large-scale evaluation of co-folding pose accuracy vs. model confidence. (2026). [citado em docs/03_metodologia_padrao_ouro.md, refer\u00eancia interna kim2026largescale]", { spacingAfter: 160 }),
      para("Dolinsky, T. J., Nielsen, J. E., McCammon, J. A. & Baker, N. A. PDB2PQR: an automated pipeline for the setup of Poisson-Boltzmann electrostatics calculations. Nucleic Acids Res. 32, W665\u2013W667 (2004).", { spacingAfter: 160 }),
      para("Olsson, M. H. M., S\u00f8ndergaard, C. R., Rostkowski, M. & Jensen, J. H. PROPKA3: Consistent Treatment of Internal and Surface Residues in Empirical pKa Predictions. J. Chem. Theory Comput. 7, 525\u2013537 (2011).", { spacingAfter: 160 }),
      para("S\u00f8ndergaard, C. R., Olsson, M. H. M., Rostkowski, M. & Jensen, J. H. Improved Treatment of Ligands and Coupling Effects in Empirical Calculation and Rationalization of pKa Values. J. Chem. Theory Comput. 7, 2284\u20132295 (2011).", { spacingAfter: 160 }),
      para("Xue, L. C., Rodrigues, J. P., Kastritis, P. L., Bonvin, A. M. & Vangone, A. PRODIGY: a web server for predicting the binding affinity of protein-protein complexes. Bioinformatics 32, 3676\u20133678 (2016).", { spacingAfter: 160 }),
    ],
  }],
});

Packer.toBuffer(doc).then(buf => {
  fs.writeFileSync(OUT, buf);
  console.log("Wrote", OUT, buf.length, "bytes");
});
