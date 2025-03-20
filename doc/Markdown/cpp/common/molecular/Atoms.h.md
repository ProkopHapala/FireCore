# Atoms.h

The `Atoms.h` header file defines the `Atoms` class, which is used to store and manipulate atomic data in molecular dynamics simulations. It provides methods for initializing, updating, and querying atomic properties.

## Includes

- `<string.h>`
- `<stdio.h>`
- `Vec3.h`
- `Mat3.h`

---

## Types (classes and structs)

---

### Class `MetaData`

`MetaData` is an abstract base class that serves as a container for additional data related to atoms. It provides a framework for storing metadata, which can be extended by derived classes.

#### properties

- `int dimensionOfDescriptorSpace`: Represents the size of the descriptor space used in machine learning or other computational tasks.

#### methods

- `void realloc(int n, bool bAtypes=true)`: Resizes the internal data structures to accommodate a new number of atoms. If `bAtypes` is true, it also resizes the array for atomic types.
- `void allocNew(int n, bool bAtypes=true)`: Allocates memory for a new set of atoms with size `n`. Similar to `realloc`, but does not check if existing data needs reallocation.
- `void dealloc(bool bAtypes=true)`: Frees the allocated memory for atomic positions and types. If `bAtypes` is true, it also frees the array for atomic types.
- `void bind(int n, int* atypes_, Vec3d* apos_)`: Binds existing data to the current object without reallocation.

---

### Class `Points`

`Points` represents a collection of points in 3D space. It is intended as a base class or utility for managing atomic positions but currently lacks implementation details and does not have any derived classes.

#### properties

- `int n`: Number of points.
- `Vec3d* ps`: Pointer to an array of `Vec3d` objects representing the coordinates of each point.

---

## Class `Atoms`

`Atoms` is a base class for systems composed of atoms, such as molecules or crystals. It manages atomic positions, types, and additional metadata.

#### properties

- `int natoms`: Number of atoms in the system.
- `int* atypes`: Array of integers representing the type of each atom.
- `Vec3d* apos`: Array of 3D vectors representing the positions of each atom. Aligned to 64 bytes for performance optimization.
- `Mat3d* lvec`: Lattice vector matrix, used in crystal systems or periodic boundary conditions (PBC).
- `double Energy`: Total energy of the system.
- `long id`: Unique identifier for the system.
- `int n0`: Number of atoms in the first part of the system (e.g., ligand).

#### methods

- `void realloc`: Resizes the internal data structures to accommodate a new number of atoms. If `bAtypes` is true, it also resizes the array for atomic types.
- `void allocNew`: Allocates memory for a new set of atoms with size `n`. Similar to `realloc`, but does not check if existing data needs reallocation.
- `void dealloc`: Frees the allocated memory for atomic positions and types. If `bAtypes` is true, it also frees the array for atomic types.
- `void bind`: Binds existing data to the current object without reallocation.
- `void copyOf`: Copies all properties from another `Atoms` instance. If memory needs to be allocated, it reallocates and copies the data.
- `int atomsFromXYZ`: Reads atomic positions and types from an XYZ format file into the current object. Optionally reallocates memory if needed.
- `void atomsToXYZ`: Writes atomic positions and types to an XYZ format file. Optionally includes the number of atoms, lattice vectors, energy, and a comment.
- `void getAABB`: Calculates the axis-aligned bounding box (AABB) for all atoms in the system.
- `Vec3d getBBcog: Returns the center of mass of the AABB.
- `void fromRigid`: Applies a rigid transformation to the atomic positions using a rotation matrix and translation vector.
- `void shift`: Shifts all atomic positions by a given vector `d`.
- `bool cellFromString`: Parses a string representing lattice vectors in XYZ format and sets the lattice vector matrix accordingly.
- `inline double measureBondLength`: Calculates the distance between two atoms with indices `ia` and `ib`.
- `inline double measureCosAngle`: Calculates the cosine of the angle formed by three atoms with indices `ic`, `ia`, and `ib`.
- `inline double measureAngle`: Calculates the actual angle between two bonds using the cosine value.

## Notes

- The class hierarchy is not fully defined in this file. `Atoms` does not inherit from any other classes.
- Memory management for arrays like `atypes`, `apos`, and `lvec` should be handled carefully to avoid memory leaks.
- The `atomsFromXYZ` function assumes the presence of a valid XYZ file with proper formatting, which may need validation in real-world applications.