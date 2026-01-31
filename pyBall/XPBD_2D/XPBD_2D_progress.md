I’ll implement the remaining features in the 2D test harness: switch to `--molecule`-only loading via [AtomicSystem](cci:2://file:///home/prokophapala/git/FireCore/pyBall/AtomicSystem.py:16:0-975:78) for `.xyz/.mol/.mol2`, add port visualization (transparent atom circles + “+” at port tips), add mouse pick+drag, and introduce `xpbd_md` vs `xpbd_relax` modes with convergence logging (momentum + bond/port error) and randomized initial distortion (pos + small complex rotation noise). Then I’ll run a `--noshow` test to verify it converges and prints diagnostics.

I’ll start by inspecting the current 2D simulator API ([upload_state](cci:1://file:///home/prokophapala/git/FireCore/pyBall/XPDB_AVBD/XPBD_2D.py:134:4-170:64), rotation buffer layout, etc.) so I can implement rotation perturbation, port visualization, and the new solver modes without guessing. Then I’ll refactor [XPBD_2D_utils.py](cci:7://file:///home/prokophapala/git/FireCore/pyBall/XPBD_2D/XPBD_2D_utils.py:0:0-0:0) to add reusable convergence diagnostics + interactive viz, and finally I’ll slim [test_XPBD_2D.py](cci:7://file:///home/prokophapala/git/FireCore/pyBall/XPBD_2D/test_XPBD_2D.py:0:0-0:0) to just CLI glue and run a `--noshow` convergence test.








I’ll implement the missing features in a reusable way (mostly in [XPBD_2D_utils.py](cci:7://file:///home/prokophapala/git/FireCore/pyBall/XPBD_2D/XPBD_2D_utils.py:0:0-0:0)), then keep [test_XPBD_2D.py](cci:7://file:///home/prokophapala/git/FireCore/pyBall/XPBD_2D/test_XPBD_2D.py:0:0-0:0) as a thin CLI wrapper. Concretely: fix port ordering (currently a real bug), add port-error diagnostics + momentum logging, add `xpbd_md`/`xpbd_relax` modes (using the new XPBD callback hook), add port visualization (transparent atoms + “+” port tips), and add mouse pick+drag similar to [test_TiledJacobi_molecules.py](cci:7://file:///home/prokophapala/git/FireCore/pyBall/XPDB_AVBD/test_TiledJacobi_molecules.py:0:0-0:0). After that I’ll run a `--noshow` convergence test.













## Implemented changes

### 1) Correct molecule loading + CLI (`--molecule`)
- [test_XPBD_2D.py](cci:7://file:///home/prokophapala/git/FireCore/pyBall/XPBD_2D/test_XPBD_2D.py:0:0-0:0) now uses:
  - `--molecule ../../cpp/common_resources/xyz/pentacene.xyz` (default)
  - Supports `.xyz/.mol/.mol2` via [pyBall.AtomicSystem](cci:2://file:///home/prokophapala/git/FireCore/pyBall/AtomicSystem.py:16:0-975:78)
- Removed the old `--geometry/--xyz` logic.

### 2) Fixed a real topology bug: **ports were wrong**
Your earlier “wrong molecules / wrong behavior” was largely because the port offsets were being written into `port_local[i, n_ports]` while the solver expects port data indexed by **neighbor slot `k`** (because the kernel iterates `k=0..3` and uses that as the port slot).  
Now in [setup_from_xyz()](cci:1://file:///home/prokophapala/git/FireCore/pyBall/XPBD_2D/test_XPBD_2D.py:86:0-189:16) ports are stored as:
- `port_local[i, k] = 0.5*(pos[j]-pos[i])` for neighbor slot `k`.

Also: bonds were rebuilt **after** node/cap reordering (previously neighs were inconsistent with the new ordering).

### 3) Initial distortion
- Position perturbation: `--perturb`
- Rotation perturbation (complex rotation): `--perturb_rot`
  - implemented as random small angle per node → complex `(cos, sin)` via [rot_from_angles()](cci:1://file:///home/prokophapala/git/FireCore/pyBall/XPBD_2D/XPBD_2D_utils.py:194:0-196:75)

### 4) Convergence diagnostics (like in your 3D tests)
Printed every 20 iters:
- **linear momentum** `|P|`
- **angular momentum** `|L|`
- **port constraint error**: `port_max`, `port_rms`

Port error is computed as the distance between matched port tips:
`(pos[i] + R_i*port_local[i,k]) - (pos[j] + R_j*port_local[j,bk])`

### 5) Solver modes (md vs relax)
`--method` choices now:
- `force`
- `xpbd_md`  (one visualization per outer loop step)
- `xpbd_relax` (visualize every XPBD Jacobi sub-iteration)

This required a small extension to the solver core:
- [XPBD_2D.step_xpbd(..., callback=None)](cci:1://file:///home/prokophapala/git/FireCore/pyBall/XPBD_2D/XPBD_2D.py:305:4-350:29) now calls `callback(itr)` after each inner iteration.

### 6) Port visualization (+ tips) + transparent atoms
[LiveViz2D](cci:2://file:///home/prokophapala/git/FireCore/pyBall/XPBD_2D/XPBD_2D_utils.py:21:0-93:29) now shows:
- atoms as **transparent circles**
- port tips as **“+” markers**
- bonds as lines

### 7) Mouse picking + dragging
Added [attach_picker_2d()](cci:1://file:///home/prokophapala/git/FireCore/pyBall/XPBD_2D/XPBD_2D_utils.py:228:0-277:15) (matplotlib callbacks), hooked in the test script.
Dragging overwrites atom position on GPU via a new core helper:
- [XPBD_2D.set_atom_pos(ia, xy)](cci:1://file:///home/prokophapala/git/FireCore/pyBall/XPBD_2D/XPBD_2D.py:207:4-213:81)

## Test run (noshow)
I ran:

```bash
python test_XPBD_2D.py --method xpbd_md --molecule ../../cpp/common_resources/xyz/H2O.xyz --iters 50 --dt 0.05 --perturb 0.2 --perturb_rot 0.2 --noshow
```

Exit code **0**, printed momentum + port error diagnostics.

## Notes / Next fixes
- In this run H2O produced `nnode=1` (oxygen only) so port error is trivially 0. That’s expected with the current heuristic (“nodes are atoms with >1 bond”). If you want H atoms also as nodes, run with `--bAllNodes`.
- There’s a PyOpenCL warning `RepeatedKernelRetrieval` (performance only). We can cache kernels later, but it’s not correctness-critical.

## Status
- All requested features are implemented and the script runs successfully in `--noshow` mode with convergence logging.

## 2026-01-31 Updates
- Fixed animation window: keep FuncAnimation alive with blocking `plt.show`, infinite frames, and stored reference; disabled blit to ensure artists render.
- Added port link lines (atom→port blue, port→neighbor magenta) using LineCollection.
- Improved picking: per-step pinning now forces picked atom position to mouse and zeros vel/omega; added `set_atom_vel`/`set_atom_omega` helpers.
- Throttling via `--viz_every`; `--noshow` runs skip per-step downloads for speed.
- Single `AtomicSystem` instantiation via `setup_from_mol`.
