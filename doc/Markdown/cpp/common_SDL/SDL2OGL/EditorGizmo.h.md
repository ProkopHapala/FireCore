# EditorGizmo.h

This header defines the `EditorGizmo` class, a Blender-like 3D transform gizmo for interactive translation, scaling, and rotation of selected points with axis constraints and selection logic. It handles mouse/keyboard events from SDL and projects motion to the screen-aligned gizmo axes.

## Includes

- SDL/OpenGL: `<SDL2/SDL.h>`, `<SDL2/SDL_opengl.h>`
- Draw/camera: `Draw.h`, `Draw2D.h`, `Draw3D.h`, `Camera.h`
- Math/geom: `Vec2.h`, `Vec3.h`, `Vec3Utils.h`, `quaternion.h`
- Picking/ray: `raytrace.h`

---

## Free Functions

- `int pickPoint(const Vec3d& ro, const Vec3d& rd, int n, Vec3d* points, double Rmax, double* Rs=0)` — Ray-pick nearest point under radius, returns index or -1.
- `inline Vec3d ray2screen(const Vec3d& v, const Vec3d& rd)` — Project 3D delta `v` into screen-parallel plane along ray `rd`.

---

## Types

- `struct Pose3d { Vec3d pos; Mat3d rot; }` — Rigid pose (origin + axes) used by the gizmo and groups.

---

## Class `EditorGizmo`

Interactive transform tool for point sets.

### Important members

- Transform/gizmo: `float Rhandle`, `float Rgizmo`, `int oglSphere`.
- Modes: `char mPick` (v/e/p), `char mTrans` (m/r/s), `char mOrig` (g/c/l), `char mPickLogic` (+/-/~), booleans `bKeyboard`, `bDragUpdate`, `bSelectAxis`.
- Axis selection: `bool axmask[3]`, `Vec3d axStart`, `Vec3d axPos`.
- Camera/mouse: `Camera* cam`, `Vec3d rd, ro`, `Vec2f pix`, `double pointSize`.
- Geometry: `int npoint, nedge`, `Vec3d* points`, `double* pointSizes`, `Vec2i* edges`.
- Selection/grouping: `std::unordered_map<int,int> selection` (point -> group), `std::vector<Pose3d> groupPose`, `int groupBrush`.
- State: `bool dragged`, `int iDebug`, `Pose3d pose`.

### Binding

- `bindPoints(int n, Vec3d* points, double* pointSizes=(double*)1)` — Bind point positions and optional radii.
- `bindEdges(int n, Vec2i* edges)` — Bind edge list for drawing/selection aids.

### Transform application

- `applyTranslation(const Vec3d& shift)` — Apply filtered translation to all selected points.
- `applyScaling(const Vec3d& sc)` — Scale selected points around `pose.pos` in `pose.rot` frame.
- `applyRotation(const Vec3d& uaxis, double phi)` — Placeholder for rotation of selected points.

### Axis handling and projection

- `clearAxMask()` — Reset axis mask to all false.
- `bool selectAxis()` — Pick closest axis handle and set `axmask` accordingly.
- `void filterByAxis(Vec3d& v)` — Zero out components of `v` not enabled in `axmask` (in gizmo frame).
- `void projectToGizmo(Vec3d& cax)` — Project gizmo origin onto screen; used for delta computation.

### Mouse and events

- `void mouseStart()` — Begin drag: optionally select axis, set `axStart/axPos`.
- `void mouseEnd()` — End/step drag: compute delta, filter by axis, apply transform according to `mTrans`.
- `void mousePick()` — Perform selection under cursor per `mPick` and `mPickLogic`.
- `void onEvent(Vec2f pix, const SDL_Event& e)` — Main event handler; updates ray from camera, handles keys/mouse.

### Keyboard bindings (when `bKeyboard`)

- Transform mode: `t` translate, `e` scale, `r` rotate.
- Pick mode: `o` vertex, `i` edge, `p` polygon.
- Origin mode: `g` global, `h` cluster, `j` local.

---

## Notes

- Motion uses screen-plane projection along camera ray for intuitive dragging, filtered by selected axes.
- Selection map stores both point index and group id to enable grouped transforms.
- Rotation mode is stubbed and can be extended to accumulate rotations around filtered axis.
