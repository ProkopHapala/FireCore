# GUI functionality to be added

Generally it make sense to be inspired by GUI of other 3D editors and molecular editors. In particular
- [Avogadro](https://avogadro.cc/) - drawing atoms (auto adjusting hydrogens) is very nice there, the rest (e.g.orientation of molecules) is rather bad
- [Blender](https://www.blender.org/) (Gizmo, vertexes=atoms, edges=bonds)
- [Godot](https://editor.godotengine.org/releases/latest/) (or other 3D game engine editors)

# Topology modification

- Add/Remove atoms 
  - probably implemented, but not deeply tested
- delete atom selection (properly adjuste bonds, Convs etc in Builder)
  - probably implemented, but not deeply tested
- add/delete bond between atoms

# Camera

- rotate, move and zoom camera (resp. viewport) by mouse or keyboard
  - make sure that boud picking and selection (box) gizmo etc. works correctly with this
- quick buttons to orient camera along axis (e.g. like in Blender)  

# Selection

- Select individual atoms and bonds by mouse picking
- Selection of atoms
  - block selection
  - add to selection 
  - remove from selection ( drag mouse + ctrl/alt/shift or other modifier ? )
- assign given selection to group of atoms (using Builder)

# Pose manipulation

- Manipulation with selection of atoms using [Gizmo](https://docs.blender.org/manual/en/latest/editors/3dview/display/gizmo.html)
  - move molecule
  - rotate molecule
  - scaling prehaps not needed
  - I already implemented some simple gizmo in `EditorGizmo.h`
- Allow to orient molecule by 3 points (atoms) or 2 bonds
   - see `Builder::orient_atoms()` and `Mat3T::fromDirUp()`
   - 1st bond (2 points) define forward-axis
   - 2nd direction (up-vector) is orthogonalized to to it ( together they define a plane)
   - 3rd direction is computed from the other two as cross product ( normal to the plane)

# Topology based pose manipulation
- rotate chemical (phenyl, -COOH) Group around bond which attach them to the backbone
  - this requires to find [graph-bridges](https://en.wikipedia.org/wiki/Bridge_(graph_theory)) - see  `splitGraphs()` and `splitGraphs()`
    - alternatively we can use `LimitedGraph.h` for more sophisticated graph-algorithms, but not sure if the synchronization of Buildder and LimitedGraph class is worth it

# Measurement
- Basic:
  - select 2 points to measure distance between atoms (aka bond length, but no need to be connected by bond)
  - select 3 points to measure angle
  - select 4 points to measure dihedral angle

# Visualization
- Color-map
  - it is quite usefull to plot bonds in some range of distances using color-scale to visualize strain in molecules (e.g. elongated bonnds red, compressed bonds blue)
    - see function `MolGUI::makeBondColoring()`
- I have implemented following functions in `MolGUI.h` for visualization of non-covalent interactions (like hydrogen bonds), but in the mean-time I forgot what they do exactly, need to test them
    ```C++
    void  showAtomGrid( char* s, int ia, bool bDraw=true );
    Vec3d showNonBond ( char* s, Vec2i b, bool bDraw=true );
    void tryPlotNonBond();
    void plotNonBondLines();
    void plotNonBondGrid();
    void plotNonBondGridAxis();
    void relaxNonBondParticles( double dt = 0.2, double Fconv = 1e-6, int niter = 1000);
    void drawParticles();
    void drawDipoleMap();
    ```

# General info
 - somewhere in the top-bar or some panel coul be these numbers:
 - number of bonds, angles and dihedrals
 - 3x3 matrix of lattice vector - `lvec` 
    - (ideally editable, which updates bond simulator in `MolWorld_sp3`/`MMFFsp3_loc` as well as `Builder`)
 - convergence of Moleculer-Dynamics simulation ( `|Fmax| < `Fconv` )
 - Invariaants of the simulation 
    - see `Vec3d cog,vcog,fcog,tqcog;` in `MolWorld_sp3`
    - total angular and linear momentum must be conserved (zero) for molecule in vacuum