# FireCore Architecture Block Diagram (C++ Components)

This diagram provides a high-level overview of the FireCore C++ codebase architecture, based on the directory structure within the \`cpp\` folder.

\`\`\`
+-----------------+     +---------------------+     +-----------------------+     +-----------------------+     +-----------------------+     +-----------------------+
|      common     | --> |        libs         | --> |         EFF         |     |   MolecularEditor     |     |         RFF         |     |     sketches_SDL      |
| (Utilities)     |     |    (Core Libraries)  |     |  (EFF Application)  |     | (Molecular Editor App)|     |  (RFF Application)  |     |    (SDL Examples)     |
+-----------------+     +---------------------+     +-----------------------+     +-----------------------+     +-----------------------+     +-----------------------+
      ^                                 ^                                         ^
      |                                 |                                         |
      +---------------------------------+                                         |
                    |                                                             |
                    v                                                             |
      +-------------------------+                                                 |
      |   common_resources    |                                                 |
      |      (Resources)      |                                                 |
      +-------------------------+                                                 |
                                                                                  |
      ^                                 ^                                         |
      |                                 |                                         |
      +---------------------------------+                                         |
                    |                                                             |
                    v                                                             v
      +---------------------+     +---------------------+     +---------------------+
      |     common_SDL    | --> |      libs_SDL     | --> |   MolecularEditor     |
      |  (SDL Utilities)  |     |  (SDL Libraries)  |     | (MolecularEditor GUI)| 
      +---------------------+     +---------------------+     +---------------------+
                                                                      ^
                                                                      |
      ^                                 ^                             |
      |                                 |                             |
      +---------------------------------+                             |
                    |                                                 |
                    v                                                 |
      +---------------------+     +---------------------+           |
      |     common_OCL    | --> |      libs_OCL     | --> |   MolecularEditor     |
      |  (OCL Utilities)  |     |  (OCL Libraries)  |     | (MolecularEditor OCL)|
      +---------------------+     +---------------------+           +---------------------+
\`\`\`

**Component Descriptions:**

*   **\`common\` (Utilities):** This is the foundation of the C++ codebase. It provides a wide range of utility libraries and header files that are used by almost all other components. This includes:
    *   Data structures and algorithms
    *   Math utilities
    *   Input/Output operations
    *   Argument parsing
    *   Scripting language integration (Lua)
    *   Testing utilities
    *   Basic data type definitions and utilities

*   **\`libs\` (Core Libraries):**  This directory contains core C++ libraries that implement key functionalities of FireCore.  Currently, the only identified library is \`libs/quadrature_lib.cpp\` (numerical integration). It's expected that other core libraries, potentially related to molecular mechanics or force fields, would reside here. These libraries depend on the \`common\` utilities.

*   **\`apps\` (Applications):** This directory houses standalone C++ applications built using the FireCore libraries. Examples include \`EFF\`, \`MolecularEditor\`, and \`RFF\`. 
    *   **\`EFF\` (Effective Force Field Application):** Likely focuses on Effective Force Field calculations, possibly for molecular simulations or property predictions.
    *   **\`MolecularEditor\` (Molecular Editor Application):** A comprehensive molecular editor with GUI, visualization, MMFF, QM/MM, and scripting capabilities. This seems to be the flagship application. Parts of Molecular Editor also depend on `common_SDL` and `common_OCL` for GUI and OpenCL functionalities respectively.
    *   **\`RFF\` (Reactive Force Field Application):** Likely focuses on Reactive Force Field (RFF) calculations, which are more advanced force fields that can model chemical reactions.
*   **Resource Dependency:**  \`common_resources\` are used as data sources by various components, especially libraries and applications.

This block diagram and description provide a preliminary architectural overview of the FireCore C++ codebase based on the observed directory structure. Further exploration of the code within each component is needed for a more detailed understanding.