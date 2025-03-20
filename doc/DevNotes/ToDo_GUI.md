# GUI Functionality to be Added

Generally, it makes sense to draw inspiration from the GUIs of other 3D editors and molecular editors. In particular:

-   [Avogadro](https://avogadro.cc/) - Drawing atoms (with auto-adjusting hydrogens) is very well-implemented there. However, other aspects like molecule orientation are less effective.
-   [Blender](https://www.blender.org/) - Consider its Gizmo, where vertices represent atoms and edges represent bonds.
-   [Godot](https://editor.godotengine.org/releases/latest/) (or other 3D game engine editors)

# Topology Modification

-   **Add/Remove Atoms**
    -   Likely implemented ( see ), but requires thorough testing.
    -  removing single atom is quite constly operation (re-ordering all atoms, confs etc.)
       => maybe it is better to have a "mask" to do some soft-removal
-   **Delete Atom Selection**
    -  Likely implemented, but requires thorough testing.
    -  see `Builder::deleteAtoms()`
    -  Should properly adjust also `bonds`, `confs` in the `Builder` class.
-   **Add/Delete Bonds** etween atoms

# Camera

-   **Camera Manipulation**
    -   Rotate, move, and zoom the camera (viewport) using the mouse or keyboard.
    -   Ensure that bounding box picking, selection, gizmos, etc., function correctly with these camera controls.
-   **Quick Camera Orientation**
    -   Implement quick buttons to orient the camera along the principal axes (similar to Blender's functionality).

# Selection

-   **Selecting individual atoms/bonds** using mouse picking.
-   **Atom Selection with boolean operations**
    -   Block selection.
    -   Add to the current selection.
    -   Remove from the current selection (using mouse drag + Ctrl/Alt/Shift or other modifier keys).
-   **Group Assignment**
    -   Assign a given selection of atoms to a group using the Builder.

# Pose Manipulation

-   **Gizmo-Based Manipulation:**
    -   Manipulate selected atoms using a [Gizmo](https://docs.blender.org/manual/en/latest/editors/3dview/display/gizmo.html).
    -   **Operations:**
        -   Translate (move) the selection.
        -   Rotate the selection.
        -   Scaling might not be necessary.
    -   A simple gizmo has already been partially implemented in `EditorGizmo.h`.
-   **Orientation by Points/Bonds:**
    -   Allow orienting a molecule using three points (atoms) or two bonds.
    -   Refer to `Builder::orient_atoms()` and `Mat3T::fromDirUp()`.
    -   **Logic:**
        -   The first bond (two points) defines the forward axis.
        -   The second direction (up-vector) is orthogonalized to the forward axis. These two define a plane.
        -   The third direction is computed as the cross product of the other two (normal to the plane).

# Topology-Based Pose Manipulation

-   **Group Rotation:**
    -   Rotate a chemical group (e.g., phenyl, -COOH) around the bond that attaches it to the backbone.
    -   This requires identifying [graph bridges](https://en.wikipedia.org/wiki/Bridge_(graph_theory)) - see `Builder::splitGraphs()` and `splitByBond()`
        -   Alternatively, `LimitedGraph.h` could be used for more sophisticated graph algorithms. However, it's uncertain if the synchronization of the Builder and LimitedGraph classes is worthwhile.

# Measurement

-   **Basic Measurements:**
    -   Select two points to measure the distance between atoms (bond length, even if not bonded).
    -   Select three points to measure the angle.
    -   Select four points to measure the dihedral angle.

# Visualization

-   **Color-Mapping:**
    -   It is useful to plot bonds within a range of distances using a color scale to visualize strain in molecules (e.g., elongated bonds red, compressed bonds blue).
    -   See the function `MolGUI::makeBondColoring()`.
-   **Non-Covalent Interaction Visualization:**
    -  I implemented following functions in `MolGUI.h` for visualizing non-covalent interactions (like hydrogen bonds), but in the meantime I forgot what they excatly do, so it needs testing:

    ```C++
    void showAtomGrid( char* s, int ia, bool bDraw=true );
    Vec3d showNonBond ( char* s, Vec2i b, bool bDraw=true );
    void tryPlotNonBond();
    void plotNonBondLines();
    void plotNonBondGrid();
    void plotNonBondGridAxis();
    void relaxNonBondParticles( double dt = 0.2, double Fconv = 1e-6, int niter = 1000);
    void drawParticles();
    void drawDipoleMap();
    ```

# General Information (GUI Display)

-   The following information could be displayed in the top bar or a separate panel:
    -   Number of bonds, angles, and dihedrals.
    -   3x3 matrix of lattice vectors (`lvec`).
        -   Ideally, it should be editable, 
          - the user updates shuld be instantly reflected in the bond simulator - i.e. synchronized with data inside `MolWorld_sp3`/`MMFFsp3_loc` and the `Builder` class instances.
    -   Convergence of the Molecular Dynamics simulation (`|Fmax|` < `Fconv`).
    -   Invariants of the simulation:
        -   See `Vec3d cog,vcog,fcog,tqcog;` in `MolWorld_sp3`.
        -   Total angular and linear momentum must be conserved (zero) 
            - at least for a molecule in a vacuum.
            - for molecule on surface it is less clear, but for molecule which is relaxed to minimum energy it should be zero as the molecule is not moving.
