#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Preparacao generica de complexo proteina-peptideo para GROMACS.

Recursos:
  - AUTO-detecta pontes dissulfeto (SG-SG distancia < 2.5 A) -> CYX
  - Configura HIS por pH: < 6.5 = HIP, 6.5-8 = HID, > 8 = HIE
  - Renumera ligante para vir apos o receptor (evita conflito pdb2gmx)
  - Atribui chain A=receptor, B=ligante
  - Numeracao serial atomica continua + TER entre cadeias
  - PRESERVA pose docada original do ligante (NAO transladar!)

Uso CLI:
  python3 prepare_complex.py \
    --receptor receptor.pdb \
    --ligante  ligante.pdb \
    --pH 7.4 \
    --out-dir prep/

  Saidas em prep/: receptor_fixed.pdb, ligante_fixed.pdb, complexo.pdb
"""

import argparse, os, sys, math


def parse_atoms(path):
    with open(path) as f:
        return [ln.rstrip("\n") for ln in f if ln.startswith(("ATOM", "HETATM"))]


def detect_disulfides(records, dist_cutoff=2.5):
    """Detecta pontes dissulfeto por distancia SG-SG."""
    cys_sg = []
    for line in records:
        if line[12:16].strip() == "SG" and line[17:20].strip() in ("CYS", "CYX"):
            resnum = int(line[22:26])
            x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
            cys_sg.append((resnum, x, y, z))
    bonded = set()
    for i, (r1, x1, y1, z1) in enumerate(cys_sg):
        for r2, x2, y2, z2 in cys_sg[i+1:]:
            d = math.sqrt((x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2)
            if d < dist_cutoff:
                bonded.add(r1); bonded.add(r2)
    return bonded


def his_form_for_pH(pH):
    """Forma de protonacao da His para um dado pH (pKa His ~ 6)."""
    if pH < 6.5:  return "HIP"   # protonada (ambos N protonados)
    if pH < 8.0:  return "HID"   # neutra, H em N-delta
    return "HIE"                 # neutra, H em N-epsilon (acima do pKa, mais comum)


def preprocess_receptor(records, cys_disulfide, his_form):
    """CYX nas pontes dissulfeto + HIS -> forma escolhida."""
    out = []
    for line in records:
        rn = int(line[22:26])
        rname = line[17:20].strip()
        aname = line[12:16].strip()
        if rn in cys_disulfide and aname == "HG":
            continue   # remove HG das CYS dissulfeto
        if rname == "CYS" and rn in cys_disulfide:
            line = line[:17] + "CYX" + line[20:]
        if rname == "HIS":
            line = line[:17] + his_form + line[20:]
        out.append(line)
    return out


def rewrite_atom(line, serial, chain, resseq):
    """Reescreve ATOM em colunas fixas, preservando coordenadas."""
    record  = line[0:6]
    aname   = line[12:16]
    altLoc  = line[16:17]
    rname   = line[17:20]
    x = line[30:38]; y = line[38:46]; z = line[46:54]
    occ = line[54:60] if len(line) >= 60 else "  1.00"
    bf  = line[60:66] if len(line) >= 66 else "  0.00"
    el  = line[76:78] if len(line) >= 78 else "  "
    if occ.strip() in ("", "0.00"): occ = "  1.00"
    if bf.strip() == "": bf = "  0.00"
    return ("{:<6s}{:>5d} {:<4s}{:<1s}{:<3s} {:<1s}{:>4d}    "
            "{:>8s}{:>8s}{:>8s}{:>6s}{:>6s}          {:>2s}").format(
        record, serial, aname, altLoc, rname, chain, resseq, x, y, z, occ, bf, el)


def reorder(records, chain, start, offset=0):
    out, s = [], start
    for line in records:
        rs = int(line[22:26]) + offset
        out.append(rewrite_atom(line, s, chain, rs))
        s += 1
    return out, s


def write_pdb(path, records, headers=None):
    with open(path, "w") as f:
        if headers:
            for h in headers: f.write(h + "\n")
        for r in records: f.write(r + "\n")
        f.write("TER\nEND\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--receptor", required=True, help="PDB do receptor")
    ap.add_argument("--ligante",  required=True, help="PDB do ligante (peptidio)")
    ap.add_argument("--pH",       type=float, default=7.4, help="pH (default 7.4)")
    ap.add_argument("--out-dir",  default="prep", help="diretorio de saida")
    ap.add_argument("--ss-cutoff", type=float, default=2.5,
                    help="cutoff distancia SG-SG (A) para detectar S-S (default 2.5)")
    args = ap.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)

    # Receptor
    print(f"[1] Lendo receptor: {args.receptor}")
    rec_raw = parse_atoms(args.receptor)
    if not rec_raw: sys.exit("ERRO: receptor vazio.")

    print(f"[2] Detectando pontes dissulfeto (cutoff {args.ss_cutoff} A)...")
    cys_ss = detect_disulfides(rec_raw, args.ss_cutoff)
    print(f"    CYS em S-S detectadas: {sorted(cys_ss)}")

    his_form = his_form_for_pH(args.pH)
    print(f"[3] Forma da His para pH {args.pH}: {his_form}")

    rec = preprocess_receptor(rec_raw, cys_ss, his_form)

    # Ligante
    print(f"[4] Lendo ligante: {args.ligante}  (pose docada preservada)")
    lig = parse_atoms(args.ligante)
    if not lig: sys.exit("ERRO: ligante vazio.")

    rec_first = int(rec[0][22:26]); rec_last = int(rec[-1][22:26])
    lig_first = int(lig[0][22:26]); lig_last = int(lig[-1][22:26])
    offset = (rec_last + 1) - lig_first
    print(f"    Receptor: residuos {rec_first}-{rec_last} ({len(rec)} atomos)")
    print(f"    Ligante : residuos {lig_first}-{lig_last} -> {lig_first+offset}-{lig_last+offset}")

    # Reordena com chain IDs
    rec_out, nxt = reorder(rec, "A", 1)
    lig_out, _   = reorder(lig, "B", nxt, offset=offset)

    # Salva
    out_rec = os.path.join(args.out_dir, "receptor_fixed.pdb")
    out_lig = os.path.join(args.out_dir, "ligante_fixed.pdb")
    out_cpx = os.path.join(args.out_dir, "complexo.pdb")

    write_pdb(out_rec, rec_out,
              [f"REMARK   Receptor pre-processado para pH {args.pH}",
               f"REMARK   {his_form} para HIS, CYX em pontes S-S: {sorted(cys_ss)}"])
    write_pdb(out_lig, lig_out,
              ["REMARK   Ligante peptidico - pose docada original preservada"])
    with open(out_cpx, "w") as f:
        f.write(f"REMARK   Complexo proteina-peptideo - pH {args.pH}\n")
        f.write(f"REMARK   Cadeia A: receptor {rec_first}-{rec_last}\n")
        f.write(f"REMARK   Cadeia B: ligante {lig_first+offset}-{lig_last+offset}\n")
        for r in rec_out: f.write(r + "\n")
        f.write("TER\n")
        for r in lig_out: f.write(r + "\n")
        f.write("TER\nEND\n")

    print(f"[5] Arquivos gerados em '{args.out_dir}/':")
    print(f"    receptor_fixed.pdb  ({len(rec_out)} atomos)")
    print(f"    ligante_fixed.pdb   ({len(lig_out)} atomos)")
    print(f"    complexo.pdb        ({len(rec_out)+len(lig_out)} atomos)")


if __name__ == "__main__":
    main()
