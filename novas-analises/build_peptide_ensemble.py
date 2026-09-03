#!/usr/bin/env python3
"""Gera ensemble de 3 conformeros iniciais de peptideo (helice-alfa, fita-beta
estendida, poliprolina-II) para uso como input flexivel do HADDOCK3, seguindo
a recomendacao de docs/frontier-tech-roadmap.md #4 (HADDOCK3 aceita/recomenda
ensemble multi-conformero de peptideo em vez de pose unica).

Roda dentro do PyMOL (env "structure" no servidor): pymol -cq build_peptide_ensemble.py -- <seq> <out_prefix>
"""
import sys
from pymol import cmd

seq = sys.argv[1]
out_prefix = sys.argv[2]

three = {
    'A':'ala','R':'arg','N':'asn','D':'asp','C':'cys','Q':'gln','E':'glu',
    'G':'gly','H':'his','I':'ile','L':'leu','K':'lys','M':'met','F':'phe',
    'P':'pro','S':'ser','T':'thr','W':'trp','Y':'tyr','V':'val',
}

def build(ss_code, name):
    cmd.reinitialize()
    cmd.fab(seq, name, ss=ss_code)
    cmd.alter(name, "chain='B'")
    cmd.save(f"{out_prefix}_{name}.pdb", name)

# helice-alfa (ss=1)
build(1, 'helix')
# fita-beta antiparalela estendida (ss=2)
build(2, 'strand')
# poliprolina-II: comeca "flat"/loop (ss=0) e ajusta phi/psi por residuo (canonico PPII: phi=-75, psi=145)
cmd.reinitialize()
cmd.fab(seq, 'ppii', ss=0)
cmd.alter('ppii', "chain='B'")
model = cmd.get_model('ppii and name CA')
resis = sorted({int(a.resi) for a in model.atom})
for r in resis:
    try:
        cmd.set_dihedral(f'ppii and resi {r-1} and name C', f'ppii and resi {r} and name N',
                          f'ppii and resi {r} and name CA', f'ppii and resi {r} and name C', -75.0)
        cmd.set_dihedral(f'ppii and resi {r} and name N', f'ppii and resi {r} and name CA',
                          f'ppii and resi {r} and name C', f'ppii and resi {r+1} and name N', 145.0)
    except Exception:
        pass
cmd.save(f"{out_prefix}_ppii.pdb", 'ppii')

print(f"OK: {out_prefix}_helix.pdb, {out_prefix}_strand.pdb, {out_prefix}_ppii.pdb")
