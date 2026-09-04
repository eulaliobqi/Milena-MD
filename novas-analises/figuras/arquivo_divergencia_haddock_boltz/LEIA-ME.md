# Arquivo — comparação HADDOCK3 × Boltz-2 (pose descartada)

Material da investigação de divergência de pose descrita em
`comparativo_vs_DN2954-GORE12T.md` (seção 7) e nas versões anteriores de
`artigo.md`. **Não faz mais parte do conjunto de figuras numeradas do
artigo** — a partir da decisão de 2026-09-03 (usar a pose Boltz-2 como
estrutura inicial de MD para os dois sistemas), `artigo.md` mostra só
figuras derivadas do Boltz-2 (`figuras/fig1_*` a `fig6_*`).

Mantido aqui por proveniência/transparência, não por uso corrente:

- `fig7_consensus_poses.{pdf,png,export.json}` — render original (HADDOCK
  padrão-ouro vs. Boltz-2 fronteira, 3 painéis) que expôs a divergência.
- `render_consensus_poses.py`, `make_fig7_consensus_poses.py` — scripts
  que geraram esse render (PyMOL no servidor + composição local).

Os números quantitativos da investigação (RMSD, distância ponta-a-ponta,
pLDDT) continuam reproduzíveis independentemente disso por
`novas-analises/scripts/diagnose_pose_divergence.py`.
