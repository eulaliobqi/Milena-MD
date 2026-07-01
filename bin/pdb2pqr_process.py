#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Limpa output do PDB2PQR para uso no GROMACS:
  - Remove hidrogênios adicionados pelo PDB2PQR (pdb2gmx os readiciona)
  - Mapeia nomes de resíduos PDB2PQR → GROMACS AMBER

Uso: pdb2pqr_process.py input_pdb2pqr.pdb output_gromacs.pdb
"""
import sys, os

# PDB2PQR AMBER → GROMACS AMBER (a maioria já é compatível)
RENAME = {
    'HISD': 'HID',   # His protonada em N-delta
    'HISE': 'HIE',   # His protonada em N-epsilon
    'HISH': 'HIP',   # His duplamente protonada
    'ASPH': 'ASH',   # Asp protonada
    'GLUH': 'GLH',   # Glu protonada
}


def process(infile, outfile):
    kept = 0
    with open(infile) as f, open(outfile, 'w') as out:
        for line in f:
            if line.startswith(('ATOM', 'HETATM')):
                atom = line[12:16].strip()
                # Remove H adicionados pelo PDB2PQR
                if atom.startswith('H') or (len(atom) > 1 and atom[0].isdigit() and atom[1] == 'H'):
                    continue
                resname = line[17:20].strip()
                new_name = RENAME.get(resname, resname)
                if new_name != resname:
                    line = line[:17] + f'{new_name:<3s}' + line[20:]
                    print(f"  Renomeado: {resname} → {new_name} (res {line[22:26].strip()})")
                kept += 1
            out.write(line)
    print(f"  {kept} átomos escritos em {outfile}")


if __name__ == '__main__':
    if len(sys.argv) != 3:
        sys.exit(f"Uso: {sys.argv[0]} input.pdb output.pdb")
    if not os.path.exists(sys.argv[1]):
        sys.exit(f"Arquivo não encontrado: {sys.argv[1]}")
    process(sys.argv[1], sys.argv[2])
