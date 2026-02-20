import os
import argparse
from pymol import cmd


def render_movie(xyz_path, outdir, prefix, width, height, dpi, ray, start, end, stride):
    if not os.path.isfile(xyz_path):
        raise FileNotFoundError(f"XYZ not found: {xyz_path}")
    os.makedirs(outdir, exist_ok=True)

    cmd.reinitialize()
    cmd.set("ray_opaque_background", 1)
    cmd.set("antialias", 2)
    cmd.set("ambient", 0.25)
    cmd.set("spec_reflect", 0)
    cmd.bg_color("white")

    # Load multi-frame XYZ as states
    cmd.load(xyz_path, "mol", multiplex=0)
    n_states = cmd.count_states("mol")
    if n_states == 0:
        raise RuntimeError("Loaded zero states from XYZ")

    start = max(1, start)
    end = n_states if end <= 0 else min(end, n_states)

    for state in range(start, end + 1, stride):
        cmd.frame(state)
        cmd.hide("everything", "mol")
        cmd.show("sticks", "mol")
        cmd.show("spheres", "elem O")  # example highlight
        cmd.color("carbon", "elem C")
        cmd.zoom("mol", buffer=2.0)
        if ray:
            cmd.ray(width, height)
        outfile = os.path.join(outdir, f"{prefix}_{state:04d}.png")
        cmd.png(outfile, width, height, dpi=dpi, quiet=1)
        print(f"wrote {outfile}")

    cmd.quit()


def main():
    ap = argparse.ArgumentParser(description="Render multi-frame XYZ to PNGs using PyMOL (headless)")
    ap.add_argument("xyz", help="Path to multi-frame XYZ movie")
    ap.add_argument("--outdir", default="pymol_frames", help="Output directory")
    ap.add_argument("--prefix", default="frame", help="PNG filename prefix")
    ap.add_argument("--width", type=int, default=800, help="Image width")
    ap.add_argument("--height", type=int, default=600, help="Image height")
    ap.add_argument("--dpi", type=int, default=200, help="PNG DPI")
    ap.add_argument("--ray", action="store_true", help="Enable ray tracing (slower, nicer)")
    ap.add_argument("--start", type=int, default=1, help="First state index (1-based)")
    ap.add_argument("--end", type=int, default=0, help="Last state index (0 means all)")
    ap.add_argument("--stride", type=int, default=1, help="Frame stride")
    args = ap.parse_args()

    render_movie(
        xyz_path=args.xyz,
        outdir=args.outdir,
        prefix=args.prefix,
        width=args.width,
        height=args.height,
        dpi=args.dpi,
        ray=args.ray,
        start=args.start,
        end=args.end,
        stride=args.stride,
    )


if __name__ == "__main__":
    main()
