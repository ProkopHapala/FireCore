# Polymer & Molecule Assembly Tools

This directory contains a suite of Python scripts for building periodic polymer chains (specifically alkanes/polyethylene), attaching functional molecules to them via anchor atoms, and combining these systems.

## Workflow Overview

The assembly process consists of three main sequential steps:
1. **Generate the Polymer Backbone:** Create a periodic zig-zag alkane chain with a specific anchor atom.
2. **Attach Molecules:** Graft side-chain molecules onto the polymer backbone at the designated anchor sites.
3. **Combine Systems:** Merge two polymer systems together, facing each other, into a single periodic cell.

---

### Step 1: Generating the Alkane Chain (`generate_alkane.py`)
Generates a periodic all-trans (zig-zag) alkane chain. It automatically computes the correct periodic boundary conditions (LVS lattice vectors) in the Y-axis. You can substitute one specific Carbon atom with a Silicon (`Si`) atom, which serves as an anchor for the next step.

**Example Usage:**
```bash
python3 ../xyz/generate_alkane.py -n 20 --lx 15.0 --lz 15.0 --si-index 10 -o my_anchor_polymer.xyz
```
*(Note: `generate_alkane.py` is located in the `../xyz/` directory).*

### Step 2: Attaching Molecules (`attach_molecules.py`)
Takes the generated polymer and a directory containing your molecules. Both the polymer and the molecules must contain exactly one Silicon (`Si`) atom acting as the anchor. The script replaces the `Si` anchors with a new `C-C` bond, calculates the missing hydrogens to maintain tetravalence, and properly aligns the molecule.

**Example Usage:**
```bash
# First, place your molecule XYZ files with 'Si' anchors into a 'molecules' directory
mkdir -p molecules

# Then run the script to attach them
python3 attach_molecules.py -p my_anchor_polymer.xyz -m molecules -o polymer_and_molecule
```

### Step 3: Combining Facing Systems (`combine_facing_systems.py`)
Takes two generated polymer-molecule systems and merges them into a single cell. The second system is automatically rotated by 180° so the attached molecules face each other. You can adjust the distance between the polymer backbones and shift them vertically/horizontally to avoid collisions or promote interactions (e.g. hydrogen bonding or intercalation).

**Example Usage:**
```bash
python3 combine_facing_systems.py \
  -s1 polymer_and_molecule/attached_MOL1.xyz \
  -s2 polymer_and_molecule/attached_MOL2.xyz \
  -dx 12.0 -dy 2.5 \
  -o combined_systems.xyz
```

**Key parameters:**
* `-dx`: Distance between the two polymer backbones (in Angstroms). Increase this for longer side-chains.
* `-dy`: Vertical shift along the polymer chain (useful for interleaving side-chains like teeth on a zipper).
* `-dz`: Horizontal shift perpendicular to the chains.