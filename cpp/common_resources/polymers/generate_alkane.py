#!/usr/bin/env python3
import argparse
import math
import os

def generate_alkane(n_carbons, lx, lz, replace_si_index=-1):
    # Počet uhlíků musí být sudý pro zachování periodicity zig-zag řetězce
    if n_carbons % 2 != 0:
        print("Varování: Počet uhlíků by měl být sudý pro správnou periodicitu. Přidávám 1.")
        n_carbons += 1

    # Parametry geometrie
    d_cc = 1.54  # délka vazby C-C v Angströmech
    d_ch = 1.09  # délka vazby C-H v Angströmech
    angle_deg = 109.47  # ideální tetraedrální úhel
    angle_rad = math.radians(angle_deg)
    half_angle = angle_rad / 2.0

    sin_half = math.sin(half_angle)
    cos_half = math.cos(half_angle)

    # Posun mezi dvěma sousedními uhlíky v ose Y a X
    delta_y = d_cc * sin_half
    delta_x = d_cc * cos_half

    atoms = []

    for i in range(n_carbons):
        # Střídání pozic pro zig-zag strukturu v ose X
        sign = 1 if i % 2 == 0 else -1
        
        c_x = sign * (delta_x / 2.0)
        c_y = i * delta_y
        c_z = 0.0

        # Identifikace atomu - možnost nahrazení za Si jako kotvu
        element = "Si" if i == replace_si_index else "C"
        atoms.append((element, c_x, c_y, c_z))

        # Pozice vodíků (míří kolmo na směr řetězce do prostoru v ose X a Z)
        h_x = c_x + sign * (d_ch * cos_half)
        h_y = c_y
        h_z = d_ch * sin_half

        atoms.append(("H", h_x, h_y, h_z))
        atoms.append(("H", h_x, h_y, -h_z))

    # Výpočet přesného rozměru periodické buňky v ose Y pro spojité navázání
    ly = n_carbons * delta_y
    
    # Posuneme atomy doprostřed definované buňky (Lx, Lz)
    for i, (el, x, y, z) in enumerate(atoms):
        atoms[i] = (el, x + lx/2.0, y, z + lz/2.0)

    # Formátování comment řádku s 'lvs' pro FireCore
    lvs_line = f"lvs {lx:.5f} 0.0 0.0    0.0 {ly:.5f} 0.0    0.0 0.0 {lz:.5f}"

    return atoms, lvs_line

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generátor periodického alkanového řetězce pro FireCore.")
    parser.add_argument("-n", "--n-carbons", type=int, default=10, help="Počet uhlíků (délka polymeru).")
    parser.add_argument("--lx", type=float, default=15.0, help="Velikost buňky ve směru X (kolmo na řetězec).")
    parser.add_argument("--lz", type=float, default=15.0, help="Velikost buňky ve směru Z (kolmo na řetězec).")
    parser.add_argument("--si-index", type=int, default=-1, help="Index uhlíku, který se má nahradit za křemík.")
    parser.add_argument("-o", "--output", type=str, default="alkane_chain.xyz", help="Název výstupního souboru.")
    
    args = parser.parse_args()

    atoms, lvs_line = generate_alkane(args.n_carbons, args.lx, args.lz, args.si_index)
    
    os.makedirs(os.path.dirname(os.path.abspath(args.output)) or ".", exist_ok=True)
    with open(args.output, 'w') as f:
        f.write(f"{len(atoms)}\n{lvs_line}\n")
        for el, x, y, z in atoms:
            f.write(f"{el:<4} {x:10.5f} {y:10.5f} {z:10.5f}\n")
            
    print(f"✅ Polymer uložen do: {args.output}")
    print(f"📐 Mřížka (lvs): {lvs_line}")