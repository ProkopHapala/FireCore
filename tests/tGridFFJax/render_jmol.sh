#!/bin/bash
# Render CHONH2/Ag(111) site structures using Jmol headless.
# Usage: bash render_jmol.sh
# Produces: top_site.png, bridge_site.png, hollow_site.png

set -e
cd "$(dirname "$0")"

for SITE in top bridge hollow; do
    echo "Rendering ${SITE}..."

    cat > _jmol_${SITE}.spt << JMOLEOF
# Load structure
load ${SITE}_scene.xyz

# Setup
set backgroundColor white
set antialiasDisplay true
set antialiasTranslucent true
select all
hide all

# Ag slab: small silver spheres
select Ag
spacefill 0.75
color [xB0B0C0]
show ball

# Molecule: larger, colored atoms with bonds
select C
spacefill 0.48
color [x202020]
show ball
show sticks

select O
spacefill 0.46
color [xDD3333]
show ball
show sticks

select N
spacefill 0.47
color [x4444EE]
show ball
show sticks

select H
spacefill 0.30
color [xFFFFFF]
show ball
show sticks

# Bond styling
set stickRadius 0.18
set bondRadiusMilliAngstroms 66

# Camera: side view (x horizontal, z up, looking along y)
moveto 1.0 {1 0 0 0} 100
rotate x -90
select all
zoom 0

# Rendering: phong lighting for depth
set lightModel phong
set ambientPercent 45
set diffusePercent 60
set specularPercent 35
set specular 0.4
set shininess 50

# Output
write image 1200 700 PNG ${SITE}_site.png
JMOLEOF

    # Run Jmol headless
    jmol -n -s _jmol_${SITE}.spt 2>/dev/null || xvfb-run -a jmol -n -s _jmol_${SITE}.spt 2>/dev/null
    rm -f _jmol_${SITE}.spt

    if [ -f "${SITE}_site.png" ]; then
        echo "  ✓ ${SITE}_site.png"
    else
        echo "  ✗ Failed to render ${SITE}"
    fi
done

echo ""
echo "Done! Now run:"
echo "  python3 plot_publication_panel.py --structure-dir ."
