#!/usr/bin/env python3
"""Render CHONH2/Ag(111) site structures using PyMOL.

Usage:  pymol -cq /tmp/ag_raw_r6_final/presentation_data/render_pymol.py
Then:   python3 plot_publication_panel.py --structure-dir /tmp/ag_raw_r6_final/presentation_data
"""
import sys
from pathlib import Path

try:
    from pymol import cmd
except ImportError:
    print("Run via: pymol -cq render_pymol.py"); sys.exit(1)

OUT_DIR = Path("/tmp/ag_raw_r6_final/presentation_data")
SITES = ["top", "bridge", "hollow"]
WIDTH, HEIGHT = 1200, 700


def render_site(site: str):
    cmd.reinitialize()
    cmd.load(str(OUT_DIR / f"{site}_scene.xyz"), "scene")

    cmd.hide("everything")

    # ---- Slab: spheres ----
    cmd.select("slab", "elem Ag")
    cmd.show("spheres", "slab")
    cmd.set("sphere_scale", 0.65, "slab")
    cmd.color("gray70", "slab")

    # ---- Molecule: ball-and-stick, enlarged ----
    cmd.select("mol", "not elem Ag")
    cmd.show("spheres", "mol")
    cmd.show("sticks", "mol")

    # Sphere sizes (molecule)
    cmd.set("sphere_scale", 0.42, "mol and elem C")
    cmd.set("sphere_scale", 0.40, "mol and elem O")
    cmd.set("sphere_scale", 0.41, "mol and elem N")
    cmd.set("sphere_scale", 0.26, "mol and elem H")

    # Colours
    cmd.color("gray10",  "mol and elem C")
    cmd.color("tv_red",  "mol and elem O")
    cmd.color("tv_blue", "mol and elem N")
    cmd.color("white",   "mol and elem H")

    # Sticks
    cmd.set("stick_radius", 0.14, "mol")
    cmd.set("stick_color", "gray30", "mol")
    cmd.set("valence", 0)

    # ---- Camera: side view ----
    # Default PyMOL: screen x=model x, screen y=model y, looking down -z
    # We want: screen x=model x, screen y=model z (up), looking along -y
    # => rotate 90° around x-axis
    cmd.reset()
    cmd.orient("all")
    cmd.turn("x", -90)       # z-axis now points up on screen
    cmd.zoom("all", 10)      # generous buffer to see everything

    # ---- Render settings (publication quality) ----
    cmd.bg_color("white")
    cmd.set("ray_trace_mode", 1)      # outline mode
    cmd.set("ray_shadow", 0)
    cmd.set("ray_trace_gain", 0.003)
    cmd.set("ray_trace_color", "black")
    cmd.set("antialias", 2)
    cmd.set("ambient", 0.50)
    cmd.set("direct", 0.35)
    cmd.set("specular", 0.25)
    cmd.set("reflect", 0.20)
    cmd.set("orthoscopic", 1)
    cmd.set("depth_cue", 1)
    cmd.set("fog_start", 0.35)
    cmd.set("ray_opaque_background", 0)  # transparent bg

    out_path = str(OUT_DIR / f"{site}_site.png")
    cmd.ray(WIDTH, HEIGHT)
    cmd.png(out_path, dpi=300)
    print(f"  {site}_site.png")


if __name__ == "__main__" or "pymol" in sys.modules:
    print("Rendering with PyMOL...")
    for s in SITES:
        render_site(s)
    print("Done!")
    cmd.quit()
