# Particle collisions by linear solvers

[Position based dynamics](https://matthias-research.github.io/pages/publications/posBasedDyn.pdf) and [projective dynamics](https://www.projectivedynamics.org/Projective_Dynamics/index.html) are very efficient methods for solving hard-degrees of freedom in molecules  like bonds ( i.e. stiff ineraction between particles with fixed topology.

Advantage of these schemes is that they are unconditionaly stable, which allows to significantly increase time-step of molecule (or other dynamical) simulations, since the hard dregrees of freedom are treated by implicit solver, rather than explicitly propagated by dynamics throu force. In force-based dynamics simulation the time step, which is generally inversly proportional to highest vibration frequency $f \approx \sqrt{k/m}$ is no loger limited by stiffnes of these hard degrees for freedom. The limitation of these methods is that linear solver can efficiently solve only linear constrains. Therefore additional non-linear (anharmonic) forces, or non-conservative forces (e.g. air drag) must be added externlly. These forces than determine the required time step.  

Here we try to extend the method to collisions between the particles (i.e. stiff short range interactions with rapidly channging topology).

## Problems:

1. The interactions are not-linear (or harmonic) unlike bonds, but are anharmonic and dissociative (i.e. like Leanrd-Jones). This means we can use linear solver only approximatively arround the minimum.

2. Due to rapid changes of the interaction pairs, it does not make much sense to explicitly construct matrixes (like in Choolesky factorization) as these matrixes would rapidly change with rapidly changing interaction topology. Therefore it makes much more sense to use iterative methods, like Jacobi or Gauss Seidel.

## Split Quasi Linear Poential

To efficiently use linear solver like Jacobi or Gauss-Seidel we will split our anharmonic potential (like Lenard-Jones) to two parts.

1. Linea (harmonic) part
2. non-linear dissociative part

Within the framework of projective dynamics we can than implement the harmonic part in the linear solver and the non-linear part as external force. The importaint assumption is that the non-linear part is small and soft (much smaller stiffness) than the harmonic part.

## The Linear poential cutoff

The linear solver tries to find such position of each particle which minimize the potential energy form all linear (harmonic) constrains. In order to function properly the linear constrains must be defined in some reasonable range around the minimum (where force from this constrain is zero). although for collision is most importaint just the repulsive part, for numerical stability and convenience we typically extend the linear region beyond the minimum to some $R_{cut} > R_{min}$. Therefore the particles-particle force can be expressed as $F=k(r_{ij}-R_{min})$ with constant stiffness $k$ along the whole interval of distance $r<R_{cut}$.

This potential is implemented into the linear solver (Jacobi or Gauss seidel).
 - At every iteration of the algorithm we check which constrains are within the range $r<R_{cut}$. Those that are we include into the solver trying to place them within distance $R_{min}$. That means if they are closer they are pushed away, if they are further, they are pulled closer.
 - If the particles are beyond $r<R_{cut}$ linear solver ignores them.

## The smoothening potential

The linear poential has serious problem: I is swithing abruptly between zero force for $r>R_{cut}$ and non-zero force $f=-k(R_{cut}-R_{min})$ when the particles cross $R_cut$ boundary. **This will lead to unstable simulations**, dificult convergence and non-physical behaviour.

For this reason we need to implement additional correction poential which transition smoothly between the harmonic one, and zero at distance. In particular we define anther cutof $R_{cut2}>R_{cut}$ for distance where the particles should not iteract by any means.

By definition this transition smoothening poential should be resposnible for disociative anharmonic nature of the non-covalent interactions, and therefore by definition must be anharmonic and concave. For this reason we cannot integrate this potential easily into the linear constrain solver. And we need to consider it as soft external poential applied to the dynamics before the constrains are solved. For this reason  this poential should be much weaked and softer than the linear part of the collision potential.

There are many functional forms how to define such correction poential. For computational performance we can limit oursefl to ponlynomical and reciprocal functions.

In following python code we implement (as an example) the simplist formula for smoothening, which is actually linear (resp. harmonic) just concave rather than convex (i.e. parabola is inverted, with negative stiffnes). Such parabola is matched to the original linar poential at $R_{cut}$ to eliminate any discontinuity in force and then original linear poential is shifted by constant $E_0$ to match the energy.  

- To ensure the force continuity at $R_{cut}$ we simply use the relation $f=-k(R_{cut}-R_{min}) = k_2(R_{cut2}-R_{cut})$
- The energy shift $E_0$ is than defined as $E_0=(k_2/2) (R_{cut2}-R_{cut})^2$


```C++
inline double getSR_x2_smooth(const Vec3d& dp, Vec3d& f, double E_min, double R_min, double R_cut, double R_cut2) {
    double r2 = dp.norm2();
    if( r2 > (R_cut2*R_cut2) ){
        f = Vec3dZero;
        return 0.0;
    }
    double r  = sqrt(r2);
    /// TODO: we should precompute these constants  k1,k2
    double d1 = R_cut  - R_min;  // harmonic region: [R_min, R_cut]
    double d2 = R_cut2 - R_cut;  // smoothing region: [R_cut, R_cut2]  
    double k1 = -2.0 * E_min / ( d1 * (d1 + d2) );
    double k2 = - k1 * (d1 / d2);
    double E, F;                 // potential and scalar force
    if (r < R_cut) {
        double x = r - R_min;
        E = 0.5 * k1 * x * x + E_min;
        F =     - k1 * x;
    } else {
        double x  = r - R_cut2;
        E = 0.5 * k2 * x * x;
        F =     - k2 * x;
    }
    f.set_mul(dp, F/r );
    return E;
}
```



---


Look inside ToDo_FastCollision.md

We want to implement the method described there. Please make a detailed plan.

To do that use

function getSR_x2_smooth() used in Forces.h
Jacobi solver and ProjectiveDynamics_d::run_LinSolve() in ProjectiveDynamics_d.cpp and ProjectiveDynamics_d.h
the whole thinkg should be implemented inside MolWorld_sp3.h in particular in MDloop
this means we should make and alterantive to run_no_omp and run_omp, something like run_projective_dynamics
The main chalange is to correctly update the neighborlist used in projective dynamics. This neighbor list should always contain all bond neighbors (see neighs inside MMFFsp3_loc ffl see MMFFsp3_loc.h and NBFF.h ) but dynamically we should also add / or remove the non-bonede atoms which are within colision distance R_cut We should not modify NBFF::neighs, we should always use them to form new neighs in ProjectiveDynamics_d before we start the Jacobi linear solver.

We should use the non-convex smothening part of the SR_x2 potential (the part betwee R_cut < r < R_cut2 ) as the external forcefield in the projective dynamics ( considering that projective dynamics is splitting scheme decomposing total interaction to the linear part solved by linar constrain solver, and non-linear part solved by external forces)

Let me know if you have any questions

Plese write comprehensive plan what steps we need to do to implement this, and note very specifically all files, classes, and variables relevant to each sub-task

---

Detailed Plan for Implementing Fast Particle Collisions using Projective Dynamics
This plan outlines the steps required to implement the fast particle collision method described in ToDo_FastCollision.md, utilizing getSR_x2_smooth(), the Jacobi solver, and ProjectiveDynamics_d. The implementation will be integrated into MolWorld_sp3.h's MDloop, creating a new run_projective_dynamics function.

I. Data Structures and Initialization (MolWorld_sp3.h):

Add necessary member variables to MolWorld_sp3:

double R_min: Minimum distance in the short-range potential.
double R_cut: Inner cutoff distance for the linear part of the short-range potential.
double R_cut2: Outer cutoff distance for the smooth transition of the short-range potential.
double E_min: Minimum energy of the short-range potential (negative value).
double k1: Stiffness of the linear part of the short-range potential.
double k2: Stiffness of the smoothing parabola. This can be pre-computed using computeMatchingParams.
std::vector<EdgeVertBond> collision_bonds: Dynamically updated list of collision pairs (atom indices and relevant parameters).
Modify init() in MolWorld_sp3.h: Initialize the new member variables above with appropriate values. Use computeMatchingParams from Forces.h to calculate R_cut and k2 based on k1, E_min, R_min, and R_cut2. These parameters should be configurable (e.g., through command line arguments or a configuration file).

computeMatchingParams (Forces.h) Integration: Ensure this function is correctly included and called during initialization in MolWorld_sp3::init().

II. Neighbor List Update (MolWorld_sp3.h and MMFFsp3_loc.h):

Create a function updateCollisionNeighbors (MolWorld_sp3.h): This function will be called within the run_projective_dynamics loop before each solver iteration. It will:

Iterate over all atoms in ffl.apos (from MMFFsp3_loc).
For each atom, iterate over all other atoms (efficiently using spatial hashing if possible).
Check the distance between each pair using dp.norm().
If the distance r is less than R_cut2, add the pair to collision_bonds. Store atom indices, distance r, and the force calculated by getSR_x2_smooth(r, R_min, R_cut, R_cut2).
Integrate updateCollisionNeighbors into run_projective_dynamics: Call this function at the beginning of each iteration of the solver.

III. Projective Dynamics Solver (ProjectiveDynamics_d.cpp and .h):

Modify prepare_LinearSystem: This function will now need to build the matrix A to include collision forces. It will still deal with bond constraints, but it will use a different method to handle the collision forces. This requires carefully considering how the matrix represents both types of interactions. Consider the following:
Sparse Matrix Representation: If PDsparse is already designed to handle sparse matrices, it might be the easiest way to add collision interactions. You could augment the existing sparse matrix with entries from collision_bonds, making the sparsity pattern dynamic. This may require adjustments to how the solver is called in run_LinSolve.
Dense Matrix Representation: If using PDmat (dense matrix), additional computational effort will be needed to fill the matrix entries that correspond to collision interactions. This might be less efficient for large numbers of particles but is simpler to implement.
Modify rhs_ProjectiveDynamics: The right-hand side of the linear system should now include the external collision forces calculated in updateCollisionNeighbors and stored in collision_bonds. Add a loop through collision_bonds. For each pair, add the pre-computed force to the appropriate entries of b.
Keep the Jacobi Solver: The existing Jacobi solver (updateJacobi_lin) is suitable for this task. Ensure its update step considers the modified b vector containing both bond and collision forces.
Add/Modify run_LinSolve: Adapt it to use the updated prepare_LinearSystem and rhs_ProjectiveDynamics. Specifically, it needs to account for dynamic changes in the neighbor list for collisions.
IV. External Force Calculation (MolWorld_sp3.h):

Create a function to calculate the external forces: This function will use getSR_x2_smooth to calculate the repulsive force between colliding particles. It will take into account the distance r between the particles and the cutoff distances (R_cut, R_cut2).
Integrate this force calculation into rhs_ProjectiveDynamics: Add the calculated force from getSR_x2_smooth to the appropriate entries in the b vector.
V. New Simulation Loop (MolWorld_sp3.h):

Create run_projective_dynamics: This function will handle the main simulation loop using projective dynamics. It will combine elements of the existing run_no_omp and run_omp (depending on your parallelization choice), incorporating the steps outlined above. This will involve:
A main loop that iterates over timesteps.
Calls to updateCollisionNeighbors at the start of each solver iteration.
Calls to prepare_LinearSystem and run_LinSolve (using the Jacobi solver).
Updating atom positions and velocities based on the solver results.
Handling of boundary conditions and other relevant aspects of the simulation.
VI. Testing and Validation:

Thorough testing: Test the implementation with various scenarios involving different numbers of particles, densities, and collision parameters. Verify that collisions are handled correctly, and the simulation remains stable.
Performance evaluation: Compare the performance of the new method with the existing methods (run_no_omp, run_omp) to assess the efficiency gain.
Validation: Compare the results with a reference simulation using a different (more accurate, but possibly slower) collision-handling method to ensure correctness.
This detailed plan provides a structured approach for implementing the proposed fast collision method. Remember to carefully handle memory allocation and deallocation to avoid errors. The use of sparse matrices will be crucial for efficiency in larger systems. Profiling and optimization will be necessary to ensure good performance.

---

Revised Plan for Implementing Fast Particle Collisions using Projective Dynamics
This revised plan addresses the concerns raised regarding the original plan. The key changes involve managing the neighbor list more efficiently within ProjectiveDynamics_d, utilizing existing bounding box structures in NBFF, and simplifying the linear solver setup.

I. Data Structures and Initialization:

Modify ProjectiveDynamics_d: Add member variables to ProjectiveDynamics_d to store collision parameters:

double R_cut: Inner cutoff distance for the linear part of the short-range potential (harmonical part).
double R_cut2: Outer cutoff distance for the smooth transition of the short-range potential.
double alphaMorse: Alpha parameter for Morse potential (from Forces.h). This is needed for correct force calculation in getSR_x2_smooth.
Initialization in MolWorld_sp3::init(): Initialize R_cut and R_cut2 in MolWorld_sp3::init(). These parameters should be configurable (command line, config file). alphaMorse is obtained from gridFF.alphaMorse.

II. Neighbor List Update (MolWorld_sp3.h and NBFF.h):

Modify NBFF::evalSortRange_BBs(): Create a modified version of NBFF::evalSortRange_BBs() called NBFF::evalSortRange_BBs_PD() (in NBFF.h). This function will:

Take a pointer to ProjectiveDynamics_d::neighs as input.
Iterate through bounding boxes as before.
For each atom pair:
Calculate the distance r using dp.norm().
If r < R_cut, add the pair to ProjectiveDynamics_d::neighs (atom indices, r, and getSR_x2_smooth parameters). Store r and the force as Quat4d parameters of neighs for later access.
If R_cut < r < R_cut2, calculate the force using getSR_x2_smooth(r, R_min, R_cut, R_cut2, alphaMorse) where R_min is taken from _mixREQ(ffl.REQs[i], ffl.REQs[j]). Add this force directly to ffl.fapos[i] and ffl.fapos[j] (no need to store).
If r >= R_cut2, no interaction.
Integrate into MolWorld_sp3::MDloop: Within MolWorld_sp3::run_no_omp and MolWorld_sp3::run_omp (and potentially new run_projective_dynamics), call ffl.evalSortRange_BBs_PD() before calling pd.run_LinSolve().

III. Projective Dynamics Solver (ProjectiveDynamics_d.cpp and .h):

Modify prepare_LinearSystem() (Remove): This function is removed entirely since we are using an iterative solver that doesn't require explicit matrix construction.

Modify rhs_ProjectiveDynamics() (No changes needed): This function remains unchanged. The modified neighbor list and force calculations in evalSortRange_BBs_PD() will correctly populate the necessary data for rhs_ProjectiveDynamics().

Use Existing Iterative Solver: Utilize run_LinSolve() and the appropriate iterative solver (e.g., updateJacobi_lin, updateGaussSeidel_lin). The external forces from non-linear interactions (R_cut < r < R_cut2) are already handled by NBFF.

IV. External Force Calculation (MolWorld_sp3.h):

No separate function needed: The external force calculation is integrated into NBFF::evalSortRange_BBs_PD().
V. New Simulation Loop (Optional):

run_projective_dynamics (Optional): A separate run_projective_dynamics function is optional. You can integrate the changes into existing run_no_omp and run_omp.
VI. Testing and Validation: Follow the same testing and validation steps as outlined in the original plan.

VII. Specific File, Class, and Variable Changes:

Files: MolWorld_sp3.h, NBFF.h, ProjectiveDynamics_d.cpp, ProjectiveDynamics_d.h, Forces.h
Classes: MolWorld_sp3, ProjectiveDynamics_d, NBFF, MMFFsp3_loc
Variables: ProjectiveDynamics_d::R_cut, ProjectiveDynamics_d::R_cut2, ProjectiveDynamics_d::neighs, NBFF::BBs, NBFF::pointBBs, MMFFsp3_loc::REQs, gridFF.alphaMorse.
This revised plan provides a cleaner and more efficient approach to integrating fast collision handling into your molecular dynamics simulation. The use of bounding boxes for efficient neighbor searching and the separation of linear and non-linear force calculations simplifies the implementation significantly. Remember to thoroughly test and validate the new implementation.

---

Revised Plan for Implementing Fast Particle Collisions using Projective Dynamics
This revised plan addresses the concerns raised regarding the original plan. The key changes involve managing the neighbor list more efficiently within ProjectiveDynamics_d, utilizing existing bounding box structures in NBFF, and simplifying the linear solver setup.

I. Data Structures and Initialization:

Modify ProjectiveDynamics_d: Add member variables to ProjectiveDynamics_d to store collision parameters:

double R_cut: Inner cutoff distance for the linear part of the short-range potential (harmonical part).
double R_cut2: Outer cutoff distance for the smooth transition of the short-range potential.
double alphaMorse: Alpha parameter for Morse potential (from Forces.h). This is needed for correct force calculation in getSR_x2_smooth.
Initialization in MolWorld_sp3::init(): Initialize R_cut and R_cut2 in MolWorld_sp3::init(). These parameters should be configurable (command line, config file). alphaMorse is obtained from gridFF.alphaMorse.

II. Neighbor List Update (MolWorld_sp3.h and NBFF.h):

Modify NBFF::evalSortRange_BBs(): Create a modified version of NBFF::evalSortRange_BBs() called NBFF::evalSortRange_BBs_PD() (in NBFF.h). This function will:

Take a pointer to ProjectiveDynamics_d::neighs as input.
Iterate through bounding boxes as before.
For each atom pair:
Calculate the distance r using dp.norm().
If r < R_cut, add the pair to ProjectiveDynamics_d::neighs (atom indices, r, and getSR_x2_smooth parameters). Store r and the force as Quat4d parameters of neighs for later access.
If R_cut < r < R_cut2, calculate the force using getSR_x2_smooth(r, R_min, R_cut, R_cut2, alphaMorse) where R_min is taken from _mixREQ(ffl.REQs[i], ffl.REQs[j]). Add this force directly to ffl.fapos[i] and ffl.fapos[j] (no need to store).
If r >= R_cut2, no interaction.
Integrate into MolWorld_sp3::MDloop: Within MolWorld_sp3::run_no_omp and MolWorld_sp3::run_omp (and potentially new run_projective_dynamics), call ffl.evalSortRange_BBs_PD() before calling pd.run_LinSolve().

III. Projective Dynamics Solver (ProjectiveDynamics_d.cpp and .h):

Modify prepare_LinearSystem() (Remove): This function is removed entirely since we are using an iterative solver that doesn't require explicit matrix construction.

Modify rhs_ProjectiveDynamics() (No changes needed): This function remains unchanged. The modified neighbor list and force calculations in evalSortRange_BBs_PD() will correctly populate the necessary data for rhs_ProjectiveDynamics().

Use Existing Iterative Solver: Utilize run_LinSolve() and the appropriate iterative solver (e.g., updateJacobi_lin, updateGaussSeidel_lin). The external forces from non-linear interactions (R_cut < r < R_cut2) are already handled by NBFF.

IV. External Force Calculation (MolWorld_sp3.h):

No separate function needed: The external force calculation is integrated into NBFF::evalSortRange_BBs_PD().
V. New Simulation Loop (Optional):

run_projective_dynamics (Optional): A separate run_projective_dynamics function is optional. You can integrate the changes into existing run_no_omp and run_omp.
VI. Testing and Validation: Follow the same testing and validation steps as outlined in the original plan.

VII. Specific File, Class, and Variable Changes:

Files: MolWorld_sp3.h, NBFF.h, ProjectiveDynamics_d.cpp, ProjectiveDynamics_d.h, Forces.h
Classes: MolWorld_sp3, ProjectiveDynamics_d, NBFF, MMFFsp3_loc
Variables: ProjectiveDynamics_d::R_cut, ProjectiveDynamics_d::R_cut2, ProjectiveDynamics_d::neighs, NBFF::BBs, NBFF::pointBBs, MMFFsp3_loc::REQs, gridFF.alphaMorse.
This revised plan provides a cleaner and more efficient approach to integrating fast collision handling into your molecular dynamics simulation. The use of bounding boxes for efficient neighbor searching and the separation of linear and non-linear force calculations simplifies the implementation significantly. Remember to thoroughly test and validate the new implementation.