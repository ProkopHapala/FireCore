# Localy Orthogonalized Orbital Forcefield

## Motivation

We want to create something between reactive forcefield, electron forcefield, and approximative quantum mechanical model of chemical bonding with very limited basis-set with molecular orbitals constrained on single atom. This should provide high speed and parallelization potential, while still capturing main essence of physical origin of chemical bonding. The hope is we create highly-transferable base-model on top of which we can build more accurate reactive force-field by systematic corrections e.g. using machine learning. 

## Nature of chemical bonding

According to Hellman-Feynman theorem, chemical bonding can be understood as attraction between the nuclei and cloud of electron density between them. Electrons clouds provide glue binding otherwise repelling nuclei. Nevertheless, electrons repel each other as well. Therefore chemical bonds are balance between attraction between nuclei-to-electrons and repulsion between nuclei-to-nuclei and electron-to-electron. Quantum mechanics bring into the picture two main phenomena. 
  1. **Delocalization of electrons** due to minimization of kinetic energy. Electrons try to be smooth slowly varying functions to minimize kinetic term $T=\left<\psi | \nabla^2 | \psi\right>$. Balance between this kinetic energy and attraction of electron to nuclei leads to shape of atomic orbitals - especially to the radial part. In our simple model we take fixed set of atomic orbital approximated by Gaussian (GTO) or exponential functions (Slater orbital, STO) with appropriate spread or decay.
  2. **Pauli repulsion** due to fermionic nature of electrons. While electrostatic interactions between delocalized electrons are rather weak, the repulsion is domineted by Pauli exclusion principle enforcing that two different electrons with the same sping cannot occupy the same quantum state. This requirement leads to orthogonality of molecular orbitals, which is tipically provided by solving by eigenstates of Hamiltonian matrix in quantum chemistry, which is the main performance bottleneck. In our model we try to eliminate need for global exact orthogonalization procedures (such as eigenvalue problem or Gram-Schmidt orthogonalization) replacing them by highly local approximations focus on two aspect of Pauli-repulsion:
       1. **hybridization** - Pauli repulsion (and to lesser degree also electrostatic repulsion) between electrons on the same atom is the origin of characteristic bonding geometry of molecules like tetrahedral (sp3), trigonal-planar (sp2), linear (sp). This can be cheaply simulated by exact orthogonalization of orbitals on the same atom. 
       2. **bond-exclusion** - Pauli repulsion between electrons on different atoms. Sigma-bonds are formed by sharing electron between two atoms. Rare excepion is boron. Pi-bonds are formed by delocalization of electron in clouds above the plane of bonds between several atoms. Why we don't see multiple sigma-bonds between same pair of atoms? Because of Pauli repulsion between electrons in these bonds. Bonding orbital contains two electrons with opposite spin. This means that it repels any other electron and tries to orthogonalize with other occupied orbitals. This is the reason why double and triple bonds are formed by pi-orbitals which are perpendicular to sigma-bonds and to each other. This is also a reason why free-electron pairs can't form bonds with orbitals of other atoms (unless they are empty). This means that for formation of bonds it is **crucial to consider occupancy of orbitals**.   
      
## The basic idea

We create approximation of quantum nature of chemical bonds by considering just interaction betwee sp-hybridized orbitals between pairs of atoms. 
Key to this idea is concept of "half-bond orbitals" (HBO). There are 4 HBOs per atom, each orthogonal to the other 3 which are created by orthogonalizing the (p_x,p_y,p_z,s) basis-functions. Mathematically each HBO $i$ of atom $A$ is linear combination of basis function $\chi \in \{p_x,p_y,p_z,s\}$ $\psi_{A,j}$ as follows: 

$$\psi_{A,i}=\sum_k c_{A,i,k} \chi_{A,k}$$, where index $k \in \{ x,y,z,s\}$

## On-site interaction - exact iterative orthogonalization

Orbitals localized on atom A $\{\psi_{A,1},\psi_{A,2},\psi_{A,3},\psi_{A,4}\}$ are kept mutually exactly orhogonal at every step of the optimization. However, this orthogonalization cannot be done by (e.g. Gram-Schmidt) as it depends on order of iteration which would introduce artifical bias to rotation of the orbitals. Instead we use symmetric iterative orthogonalization methods like Jacobi rotation which are order independent.

```C++ 
template<int n>
void orthonormalize_atom(vec<n>* orbitals double tol=1e-3, int nMaxIter=10, double damp=1.0 ){
    vec<n> forces[n];
    for(int i=0;i<nMaxIter;i++){  
        double res2=0;
        // --- evaluate forces (rotations?) 
        for(int i=0;i<n;i++){
            vec<n>* oi = orbitals[i];
            for(int j=i+1;j<n;j++){
                vec<n>* oj = orbitals[j];
                double cij = dot(oi,oj);
                if(cij<tol) continue;
                forces[i] += oi*cij;
                forces[j] += oj*cij;
            }
        }
        // --- apply rotations 
        for(int i=0;i<n;i++){
            res2 += dot(forces[i],forces[i]);
            orbitals[i] -= forces[i]*damp;
            orbitals[i].normalize();
        }
        if(res2<tol) break;
    }
}
```

## Inter site interactions

In further discussion we will heavily rely on calculation of overlaps $S_{ab}=\left< \psi_a | \psi_b \right>$ between orbitals $\psi_a$ and $\psi_b$ localized on different atoms $A$ and $B$. Since the orbitals $\psi_a$ are linear combination of 4 basis function $\{\chi_x, \chi_y, \chi_z, \chi_s\}$ such overlap would comprise of 16 integrals $\left< \chi_{A,k} | \chi_{B,l} \right>$, where $k,l \in \{x,y,z,s\}$. Nevertheless if we rotate the coordinate system so that the line connecting atoms A-B $r_{AB}$ is along $z$-axis, we can simplify the overlap matrix considerably due to symmetry (odd/even functions along each cartesian axis) to only 5 elements $\{S_{xx},S_{yy},S_{zz},S_{ss},S_{sz}\}$ (assuming overlap matrix is symmetric, $S_{ab}=S_{ba}$). The overlpa matrix looks like this:
```
    x  y  z  s
   ------------
x | xx 0  0  0 
y | 0  yy 0  0 
z | 0  0  zz sz 
s | 0  0  sz ss
```    

## Electrostatics

* Electrostatic interaction betwee nuclei is solved simply by Coulomb potential $V_{nn}=\frac{Q_A Q_A}{r_{AB}}$, where charge of $A$-th atom $Q_A = Z_A - n^{CORE}_A$ is atomic number of atom $Z_A$ decreased by number of core electrons $n^{CORE}_A$.
* Attraction betwee valence electrons of other atoms and nuclei is the major component of the bonding energy, and the main reason for formation of bonds. It is evaluated according to proper quantum mechanical manner as interaction of potential of atom $V_A$ of atom $A$ and density $\rho_i=|\psi_i|^2$ of electon $i$, as follows: 

$$E_{Ai} = \int_r \rho_i(r) V_A(r) = \left< \psi_i \left| \frac{Q_A}{|r_{ij}|} \right| \psi_i \right> $$ 

* Notice that basis function $\chi_{A,k}$ are mutually orhogonal, $\left< \chi_{A,k} | \chi_{A,k'} \right> = \delta_{kk'}$, therefore the denisty $\rho_i(r) = |\sum_{ik} c_{ik}^2 \chi_{A,k}(r)|^2$ can be represented just by the trace $\rho_i(r) = \sum_{ik} c_{ik}^2 |\chi_{A,k}(r)|^2 $. The integrals  $ V_A(r,\theta) = \int_r |\chi_{A,k}(r,\theta)|^2 / |r|  $ can be evaluated easily. Nice analytical formulas exist for Gaussian basis-set, for other basis function can be evaluated efficiently using numerical integration in cylindrical coordinates.

## Electron-electron interaction

* Each half-bond orbital has floating point occupation number $f_i$ which determines strenght of its intaction with orbitals on the other atoms.
   * If we want to be rigorous, we would introduce separate occupation numbers $f_{i,\uparrow}$ and $f_{i,\downarrow}$ for spin up and down. Nevertheless, this would introduce some complexty we do not want to deal with at the moment.  
   * Instead we indroduce concept of *other-electron* occupation number $o_i \in [0.,1.]$ which track total ammount of overlap with oribitals on other atoms without distingushing their spin.
      * For hybridized orbital $\psi_a$ on atom $A$ we can evalute it as 
        $$o_a = \sum_i S_{ai} = \sum_i \left< \psi_a | \psi_i \right>$$ for all orbitals $i$ localized on other atoms.

      * We can then evaluate pentalty energy function for over-saturated bonds.
        * **orbital over-saturation** - $E^{os}_a = K_{os} (o_a - 1)^2$ 
           * This energy term enforces accepting just 1 electron from other atoms (froming single bond)
               * Notice that this energy term can be attractive if $o_a<1$. 
               * It may make sense to use separate stiffness constant $K_{os}$ for attractive $o_a<1$ and repulsive $o_a>1$ regime.
           * We may modify this depending on occupation of the orbital itself. 
              * For example for electron pairs which are already double occupied we can use $E^{os}_a = K_{os} o_a^2$ 
              * while for empty orbitals we can use $E^{os}_a = K_{os} (o_a-2)^2$. 
              * In effect we can use general formula $E^{os}_a = K_{os} (o_a+o^0_a-2)^2$, where $o^0_a$ is base occupation of the orbital $a$ itself.    
        * **bond over-saturation**    - $E^{bs}_{ab} = K_{bs} (o_a - S_{ab})(o_b - S_{ab})$  
           * Alternatively  $E^{bs}_{ab} = K_{bs} (1-S_{ab}/(o_a o_b))$
           * This term is always repulsive.
* **Electron numbering approximation** - The method described above consider overlap of all orbitals localized on other atoms. This could be rather costly - consider we have $N$ atoms in the neighborhood of the atom $A$. This means we need to consider $16N$ orbital overlaps $S_{ai}$, and for each such overlap we calculate (or interpolate from integral table) 16 overlaps between basis functions $\chi_{A,k}$ and $\chi_{B,l}$ (although this reduces when considering simmetries to just 5 interactions). This means we need to do 5-interpolations (or analytical evaluations) and 80 coefficiet multiplications per pair of atoms which is significat. 
  * We may still reduce this complexity if we consider than only orbitals with the same ordering number can form a bond. This affect calculation of the occupation number $o_a$ which is now computed only for the same slots of
To reduce this computational complexity  
* **Electrostatic interaction between electrons** is not essential for description of chemical bonding, as it can be easily captured re-parametrization of the Pauli-repulsion bonding interaction described above. But to capture long-range interactions including hydrogen-bonds it is useful to introduce electrostaic interaction betwee electrons explicitly. 

## Capping polarized hydrogen

* While normal atoms (notably C, O, N) have 4-orbitals (hybridization sp3,sp2,sp), could be described by a single s orbital. However such simple description could cause several problems:
  * The atom ca be bonded by several other atoms (as it is accesible form all directions)
  * The electron cloud should deflect toward the bond, which is not possible with single s-orbital.
     * This is importaint for formation of dipole which create hydrogen-bonding.
* To resolve these issues we sould still formed the orbital on the hydrogen as linear combination of the 4 basis functions $\{p_x,p_y,p_z,s\}$. 
        