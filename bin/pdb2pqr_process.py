#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Limpa output do PDB2PQR para uso no GROMACS:
  - Remove hidrogênios adicionados pelo PDB2PQR (pdb2gmx os readiciona)
  - Mapeia nomes de resíduos PDB2PQR (sempre AMBER -- chamado com
    `--ff AMBER --ffout AMBER` no PREPARE_PH, independente da força de
    campo final da MD) -> convenção de resíduo da força de campo alvo.

Uso: pdb2pqr_process.py input_pdb2pqr.pdb output_gromacs.pdb [ff_style]
     ff_style: amber (default) ou charmm
"""
import sys, os

# PDB2PQR (--ff AMBER) → GROMACS AMBER: a maioria já é compatível, só os
# estados de protonação não-padrão (HISD/E/H, ASPH, GLUH) precisam mapear
# pro nome de resíduo AMBER usado no .rtp do GROMACS (HID/E/P, ASH, GLH).
RENAME_AMBER = {
    'HISD': 'HID',   # His protonada em N-delta
    'HISE': 'HIE',   # His protonada em N-epsilon
    'HISH': 'HIP',   # His duplamente protonada
    'ASPH': 'ASH',   # Asp protonada
    'GLUH': 'GLH',   # Glu protonada
}

# PDB2PQR (--ff AMBER) → CHARMM36: PDB2PQR só sabe emitir nomenclatura
# AMBER-flavored (não tem --ffout CHARMM utilizável aqui). Mapeia pro nome
# CHARMM final -- cobre AMBIGUAS as duas formas que pdb2pqr pode emitir:
# verificado empiricamente (2026-09-04, receptor DN2954/charmm36) que ele
# escreve o nome AMBER FINAL direto (HID/HIE, não a forma intermediária
# HISD/HISE que o dict original assumia sem checar) -- por isso as chaves
# HID/HIE/HIP também estão aqui, não só HISD/HISE/HISH. ASPH/GLUH nunca
# apareceram num teste real (sistema sem ASP/GLU protonado em pH 8.2), mas
# ambas as formas ficam cobertas por segurança. CYX (que o PRÓPRIO pdb2pqr
# já emite pra dissulfeto em --ff AMBER, antes deste script rodar) -> CYS2,
# já que "CYX" não existe no .rtp CHARMM (confirmado em produção: grompp
# abortou com "Residue 'HIE' not found" na 1a tentativa CHARMM36 -- história
# completa no commit que introduziu --ff-style).
RENAME_CHARMM = {
    'HISD': 'HSD', 'HID': 'HSD',
    'HISE': 'HSE', 'HIE': 'HSE',
    'HISH': 'HSP', 'HIP': 'HSP',
    'ASPH': 'ASPP', 'ASH': 'ASPP',
    'GLUH': 'GLUP', 'GLH': 'GLUP',
    'CYX':  'CYS2',
}

RENAME_TABLES = {'amber': RENAME_AMBER, 'charmm': RENAME_CHARMM}


def process(infile, outfile, ff_style='amber'):
    rename = RENAME_TABLES[ff_style]
    kept = 0
    with open(infile) as f, open(outfile, 'w') as out:
        for line in f:
            if line.startswith(('ATOM', 'HETATM')):
                atom = line[12:16].strip()
                # Remove H adicionados pelo PDB2PQR
                if atom.startswith('H') or (len(atom) > 1 and atom[0].isdigit() and atom[1] == 'H'):
                    continue
                # Nome do residuo: colunas 18-21 (4 chars) -- CYS2/GLUP/ASPP
                # (CHARMM) tem 4 letras, não cabem no campo padrão de 3
                # colunas. A coluna 21 é espaço em branco no PDB estrito
                # (antes do chain ID em 22), então usar as 4 colunas aqui
                # não desloca nada pros nomes de 3 letras (mesmo raciocínio
                # de prepare_complex.py).
                resname = line[17:21].strip()
                new_name = rename.get(resname, resname)
                if new_name != resname:
                    line = line[:17] + f'{new_name:<4s}' + line[21:]
                    print(f"  Renomeado: {resname} → {new_name} (res {line[22:26].strip()})")
                kept += 1
            out.write(line)
    print(f"  {kept} átomos escritos em {outfile}")


if __name__ == '__main__':
    if len(sys.argv) not in (3, 4):
        sys.exit(f"Uso: {sys.argv[0]} input.pdb output.pdb [amber|charmm]")
    if not os.path.exists(sys.argv[1]):
        sys.exit(f"Arquivo não encontrado: {sys.argv[1]}")
    ff_style = sys.argv[3] if len(sys.argv) == 4 else 'amber'
    if ff_style not in RENAME_TABLES:
        sys.exit(f"ff_style desconhecido: {ff_style} (use amber ou charmm)")
    process(sys.argv[1], sys.argv[2], ff_style)
