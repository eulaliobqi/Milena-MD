#!/usr/bin/env python3
"""Renderiza poses consenso (HADDOCK3 padrao-ouro vs Boltz-2 fronteira) para os
dois sistemas DN2954 x GORE1-2T / GORE1-2T(GGS)3, no mesmo referencial de
camera, para figura de publicacao do Bloco 1 de novas-analises.

Roda dentro do PyMOL (env "structure" no servidor):
  pymol -cq render_consensus_poses.py -- <out_dir>

Alinha cada pose (HADDOCK, Boltz-2 model_0) pela cadeia A (receptor) contra
uma estrutura de referencia comum (data/DN2954-receptor.pdb), garantindo que
os dois paineis (peptideos diferentes, mesma tripsina DN2954) usem a mesma
orientacao de camera para comparacao visual direta.
"""
import sys
from pymol import cmd

argv = sys.argv[sys.argv.index("--") + 1:] if "--" in sys.argv else sys.argv[1:]
out_dir = argv[0] if argv else "."

BASE = "/home/eulalio/gromacs/Milena-MD/novas-analises"
REF = f"{BASE}/data/DN2954-receptor.pdb"

SYSTEMS = {
    "gore12t": {
        "haddock": f"{BASE}/DN2954-GORE1-2T/haddock/best_pose_haddock.pdb",
        "boltz": f"{BASE}/DN2954-GORE1-2T/boltz2/boltz_out/boltz_results_dn2954_gore1-2t/predictions/dn2954_gore1-2t/dn2954_gore1-2t_model_0.pdb",
    },
    "ggs3": {
        "haddock": f"{BASE}/DN2954-GORE1-2T-GGS3/haddock/best_pose_haddock.pdb",
        "boltz": f"{BASE}/DN2954-GORE1-2T-GGS3/boltz2/boltz_out/boltz_results_dn2954_gore1-2t-ggs3/predictions/dn2954_gore1-2t-ggs3/dn2954_gore1-2t-ggs3_model_0.pdb",
    },
}

TETRAD_RESI_HADDOCK = "79+126+222+216"  # numeracao PDB original, sempre usada no objeto haddock

# Okabe-Ito, normalizado 0-1
COL_HADDOCK = (0.000, 0.447, 0.698)  # azul
COL_BOLTZ = (0.835, 0.369, 0.000)    # vermelhao
COL_TETRAD = (0.000, 0.620, 0.451)   # verde-azulado
COL_RECEPTOR = (0.780, 0.780, 0.780)

cmd.set("ray_opaque_background", 1)
cmd.bg_color("white")
cmd.set("antialias", 2)
cmd.set("cartoon_fancy_helices", 1)
cmd.set("ray_trace_mode", 0)
cmd.set("specular", 0.15)
cmd.set("cartoon_transparency", 0.0)
cmd.set("depth_cue", 0)
cmd.set("ray_trace_fog", 0)
cmd.set("fog", 0)

cmd.load(REF, "ref")
cmd.hide("everything", "ref")
cmd.show("cartoon", "ref and polymer")
cmd.color("grey80", "ref")

objs = {}
for tag, paths in SYSTEMS.items():
    hobj = f"haddock_{tag}"
    bobj = f"boltz_{tag}"
    cmd.load(paths["haddock"], hobj)
    cmd.load(paths["boltz"], bobj)

    cmd.align(f"{hobj} and chain A and polymer", "ref and chain A and polymer", object=f"aln_h_{tag}")
    cmd.align(f"{bobj} and chain A and polymer", "ref and chain A and polymer", object=f"aln_b_{tag}")
    cmd.delete(f"aln_h_{tag}")
    cmd.delete(f"aln_b_{tag}")

    cmd.hide("everything", hobj)
    cmd.hide("everything", bobj)

    # peptideo HADDOCK
    cmd.show("cartoon", f"{hobj} and chain B and polymer")
    cmd.color("0x0072B2", f"{hobj} and chain B")

    # peptideo Boltz-2 (so a cadeia B; receptor do objeto boltz fica oculto, usa-se "ref")
    cmd.show("cartoon", f"{bobj} and chain B and polymer")
    cmd.color("0xD55E00", f"{bobj} and chain B")

    objs[tag] = (hobj, bobj)

# triade catalitica (numeracao PDB original) no objeto de referencia comum
tetrad_sel = "ref and chain A and resi " + TETRAD_RESI_HADDOCK
cmd.show("sticks", f"{tetrad_sel} and not name C+N+O")
cmd.color("0x009E73", tetrad_sel)
cmd.util.cnc(tetrad_sel)

# camera unica: enquadra o receptor inteiro (nunca cortado) + o trecho do
# peptideo proximo da interface (<=15 A do receptor) nos dois sistemas.
# GORE1-2T(GGS)3 e um peptideo de 75 aa que, na pose HADDOCK, sai da
# interface como helice praticamente reta e longa (~110 A) -- mostrar essa
# cauda inteira no mesmo enquadramento do painel comparativo espremeria o
# receptor a um ponto ilegivel, entao o enquadramento aqui prioriza a
# regiao de contato (achado a registrar no texto, nao escondido: ver
# consensus_ggs3_full.png para a pose completa sem esse corte de camera).
near_chainB = " or ".join(
    f"(({h} or {b}) and chain B and byres (all within 15 of (ref and chain A and polymer)))"
    for h, b in objs.values()
)
tetrad_only = tetrad_sel
cmd.orient(f"(ref and polymer) or {near_chainB}")
cmd.zoom(f"(ref and polymer) or {near_chainB} or {tetrad_only}", buffer=5)
saved_view = cmd.get_view()

for tag, (hobj, bobj) in objs.items():
    cmd.show("cartoon", "ref and polymer")
    for other_tag, (oh, ob) in objs.items():
        show = other_tag == tag
        cmd.set("cartoon_transparency", 0.0 if show else 1.0, f"{oh} and chain B")
        cmd.set("cartoon_transparency", 0.0 if show else 1.0, f"{ob} and chain B")

    cmd.set_view(saved_view)
    cmd.ray(2000, 2000)
    cmd.png(f"{out_dir}/consensus_{tag}.png", dpi=600)

# Render "full" do GGS3: mesma cena, sem restringir a camera a <=15 A do
# receptor -- documenta a pose HADDOCK inteira (helice ~75 aa quase reta),
# achado citado na legenda em vez de cortado da figura principal.
tag = "ggs3"
hobj, bobj = objs[tag]
for other_tag, (oh, ob) in objs.items():
    show = other_tag == tag
    cmd.set("cartoon_transparency", 0.0 if show else 1.0, f"{oh} and chain B")
    cmd.set("cartoon_transparency", 0.0 if show else 1.0, f"{ob} and chain B")
cmd.orient(f"(ref and polymer) or (({hobj} or {bobj}) and chain B and polymer)")
cmd.zoom(f"(ref and polymer) or (({hobj} or {bobj}) and chain B and polymer)", buffer=4)
cmd.ray(2600, 1400)
cmd.png(f"{out_dir}/consensus_ggs3_full.png", dpi=600)

print("OK")
