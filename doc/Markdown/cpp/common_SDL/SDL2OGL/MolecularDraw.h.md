# MolecularDraw.h

Utility drawing helpers for molecular visualization and field-based rendering using OpenGL. Provides color maps, neighbor/trajectory visualization, and substrate/isosurface rendering over grid force fields.

## Includes

- Data/utils: `AtomicConfiguration.h`, `molecular_utils.h`

---

## Color helpers

- `void colorRB(float f)` — Red–blue scale around gray for scalar `f`.
- `void colorRBH(float f, float h)` — Red–blue scale modulated by hue/height `h`.
- `void colorBW(float f)` — Black–white grayscale around mid-gray.

---

## Debug/printing

- `void printPoses(int n, double* poses)` — Print array of packed poses (8 doubles per item).

---

## Neighborhood and mapping

- `void drawMapedPoints(const FastAtomicMetric& D, int itest)` — Show atoms, cell coverage, and sphere around selected atom; debug mapped cell indices.
- `void drawNeighs(const FastAtomicMetric& D, Vec3d pos)` — Visualize neighbor search at `pos`: cutoff sphere, host grid cell, and connecting lines to found neighbors.

---

## Grid-force interaction and trajectories

- `void drawPPRelaxTrj(int n, double dt, double damp, GridFF& gff, Vec3d pos, Quat4f PRQ)` — Integrate a probe particle in a force field; plot relaxation trajectory.
- `void drawGridForceAlongLine(int n, GridFF& gff, Vec3d pos0, Vec3d dpos, Quat4f PRQ, double fsc)` — Sample force vectors along a line and draw arrows.

---

## Planes and substrates

- `void plotSurfPlane(Vec3d normal, double c0, Vec2d d, Vec2i n)` — Draw a reference plane grid defined by normal and extents.
- `int renderSubstrate_(const GridShape& grid, Quat4f* FF, Quat4f* FFel, double isoval, bool sign, float sclr=1.0)` — Build triangle strips of an isosurface z-height from vector field; color by electrostatic component.
- `int renderSubstrate_new(const GridFF& gff, Vec2d zrange, double isoval, Quat4d PLQ, double sclr, bool bErrNan=false)` — Robust isosurface finder between z-planes using `GridFF::findIso`; colors and normals from PL/Q contributions.
- `void viewSubstrate(int nx, int ny, int isoOgl, Vec3d a, Vec3d b, Vec3d pos0=Vec3dZero)` — Tile a compiled isosurface display list over a lattice of translations.
- `void viewSubstrate(Vec2i nxs, Vec2i nys, int isoOgl, Vec3d a, Vec3d b, Vec3d pos0=Vec3dZero)` — As above with explicit index ranges.

---

## Notes

- Immediate-mode OpenGL is used for clarity; heavy meshes should be precompiled into display lists (see `viewSubstrate`).
- `GridFF` methods supply interpolation, force accumulation, and iso-finding; these helpers focus on visualization.
- Neighbor visualization assumes `FastAtomicMetric` exposes `ruler`, `pos`, `Rcut`, and neighbor search API.
