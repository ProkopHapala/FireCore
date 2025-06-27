**Effective Core Pseudopotentials (ECP) in eFF**

The original formulation of the electron Force Field (eFF) describes an *N*-electron wave function as a Hartree product of one-electron floating spherical Gaussian (FSG) wave packets. This formulation was limited to low-Z elements, primarily those with electrons of predominant s-character. The intrinsic limitation arises from the spherical symmetry of the underlying Gaussian basis functions. For atoms with valence electrons possessing higher angular momenta, such as p-block elements, the FSG representation does not adequately capture certain interactions, particularly between the core and valence electrons, due to the absence of nodal structures.

To overcome this limitation and extend eFF to support higher-Z elements, a formal set of potential form extensions known as **Effective Core Pseudopotentials (ECPs)** were introduced. In this eFF–ECP scheme, the interaction between the core and valence electrons is replaced with a potential energy based on their overlap.

**How the ECP Model Works in eFF**

In the ECP framework, a **model represents the core electrons of an atom together with the nucleus as a single pseudo particle with wave-like behavior**. This pseudo-core particle interacts with valence electrons, other nuclei, and other pseudo-cores through effective core pseudopotentials. Nuclei are still treated as classical point charges.

The standard electrostatic interactions in eFF include nucleus–nucleus (ENN), electron–electron (Eee), and nucleus–electron (ENe) terms. The ECP scheme reformulates and parameterizes the Pauli energy term (EPR) for pseudo particles replacing the core electrons and nucleus, and adjusts the classical electrostatic energies between the pseudo-core and valence electrons (core-elec), nuclei (core-nuc), and other pseudo-core (core-core) particles.

The charge for a pseudo-core atom is given by the number of outer shell electrons. A pseudo-core should be described with a mass that includes the corresponding nuclear mass plus all the core electrons (but not outer shell electrons). Its radius is equivalent to that of a corresponding minimized full-electron system.

**Pauli Potential Modification with ECP**

A key component of the ECP extension is the reformulation of the spin-dependent Pauli repulsion potential term (EPR). In the ECP model representation, the new Pauli potentials are designed according to the relation $E_{PR} \propto S^2$, where $S$ is the overlap between two Gaussians: one representing the core and the other an interacting valence electron.

By choice, two different types of overlaps are defined in this ECP representation:
1.  **s-s overlap**: For an s-type pseudo-core interacting with an s-type valence electron.
2.  **s-p overlap**: For an s-type pseudo-core interacting with an s-type Gaussian representing one of the lobes of a p-type valence electron.

The corresponding functional forms for the ECP-modified Pauli potential are provided:

*   For core interactions with **s-type valence electrons (ECP\_s)**:
    $E_{Pauli(ECP\_s)} = p_1\exp\left(-\frac{p_2r^2}{p_3+s^2} \right)$
    Here, $r$ is the distance between the s-type pseudo-core and the interacting s-type valence electron. $s$ is the size of the valence electron. $p_1$ corresponds to the pseudo-core wave function amplitude, $p_2$ to the pseudo-core wavefunction decay factor, and $p_3$ to the square of effective pseudo-core particle size.
*   For core interactions with **p-type valence electrons (ECP\_p)**:
    $E_{Pauli(ECP\_p)} = p_1\left( \frac{2}{p_2/s+s/p_2} \right)\left( r-p_3s\right)^2\exp \left[ -\frac{p_4\left( r-p_3s \right)^2}{p_5+s^2} \right]$
    Here, $r$ corresponds to the distance between the s-type pseudo-core and the s-type Gaussian representing one of the lobes of a p-type valence electron. $s$ is the size of the valence electron. $p_1$ corresponds to the pseudo-core wave function amplitude, $p_4$ to the pseudo-core wavefunction decay factor, and $p_5$ to the square of effective pseudo-core particle size. For the s–p case, $p_3$ corresponds to an off-center measure, and $p_2$ to a second effective size that adjusts the overlap amplitude.

These parameters ($p_1, p_2, p_3$ for s-s; $p_1, p_2, p_3, p_4, p_5$ for s-p) are optimized against quantum mechanical calculations of small representative molecules.

**Elements Supported by ECP**

ECPs are supported and validated for most of the second and third row elements of the p-block. Predefined parameters are provided for **C, N, O, Al, and Si**. Al and Si can be accurately described using the s-s ECP form, while C, N, and O require the s-p form due to their more complex and dominant p-type interactions, including multiple bonds and lone pairs.

**Benefits and Limitations of ECP in eFF**

*   **Benefits**:
    *   Enables accurate description of p-block elements and systems with increasingly non-spherical electrons.
    *   Allows for the description of complex bonding structures (e.g., multiple bonds and lone pairs).
    *   Leads to a significant reduction in the total number of degrees of freedom by removing core electrons.
    *   Filters out high-frequency modes from core electrons, resulting in larger integration time steps and enabling longer simulation timescales (10s of nanoseconds for millions of atoms).
    *   Maintains the scalability capabilities of the all-electron eFF methodology.
    *   Improves geometry predictions and bonding energy calculations for elements like Al, Si, C, and O compared to the all-electron eFF.
*   **Limitations/Challenges**:
    *   Requires further extension to d-elements (e.g., using angular momentum projection operators).
    *   Needs implicit hybridization support.
    *   Requires appropriate accounting for conjugation.
    *   Development of new types of ECPs with parameters that retain physical meaning is desired.

## Particle Types and Interaction Matrix in eFF-ECP

### Spin Values and Particle Types
In the eFF-ECP implementation, particles are distinguished by their spin values, which determine their type and interaction behavior:

| Spin Value | Particle Type | Description |
|------------|---------------|-------------|
| 0 | Nucleus (N)               | Classical point charge representing an atomic nucleus |
| 1 | Electron (E)              | Valence electron represented as a Gaussian wave packet |
| 2 | Fixed Core (FC)           | Core electrons treated as a fixed potential |
| 3 | s-type Pseudo-Core (PC_s) | Pseudo-particle representing core electrons with s-symmetry |
| 4 | p-type Pseudo-Core (PC_p) | Pseudo-particle representing core electrons with p-symmetry |

## Detailed Interaction Matrix

### Abbreviations
- **ENN(i,j)** = `ElecNucNuc(q[i]*q[j], rc, &ecoul, &fpair)`
- **ENC(i,j)** = `ElecNucElec(q[i], rc, eradius[j], &ecoul, &fpair, &force)`
- **EEE(i,j)** = `ElecElecElec(rc, eradius[i], eradius[j], &ecoul, &fpair, &f1, &f2)`
- **EEC(i,j)** = `ElecElecCore(q[j], rc, eradius[j], eradius[i], &ecoul, &fpair, &f1)`
- **ECE(i,j)** = `ElecCoreElec(q[j], rc, eradius[j], eradius[i], &ecoul, &fpair, &f1)`
- **ECN(i,j)** = `ElecCoreNuc(q[j], rc, eradius[j], &ecoul, &fpair, &f1)`
- **ECC(i,j)** = `ElecCoreCore(q[i]*q[j], rc, eradius[i], eradius[j], &ecoul, &fpair)`
- **PEE0(i,j)** = `PauliElecElec(0, rc, eradius[i], eradius[j], ...)`
- **PEE1(i,j)** = `PauliElecElec(1, rc, eradius[i], eradius[j], ...)`
- **PCE(i,j)** = `PauliCoreElec(...)` (s-type)
- **PCP(i,j)** = `PauliCorePElec(...)` (p-type)

*Repeated calls → multiple electrons or spin channels.*

| i \ j | N (0)             | E (1)                                                  | FC (2)                                                                               | PCs (3)                                                 | PCp (4)                                                 |
|-------|-------------------|--------------------------------------------------------|--------------------------------------------------------------------------------------|---------------------------------------------------------|----------------------------------------------------------|
| **0** | ENN(0,0)×1       | ENC(0,1)×1                                            | ENN(0,2)×1                                                                            | ENN(0,3)×1                                              | ENN(0,4)×1                                               |
| **1** | ENC(1,0)×1       | EEE(1,1)×1 + PEE0(1,1)×1 + PEE1(1,1)×1                  | ENC(1,2)×1 + EEE(1,2)×2 + PEE0(1,2)×1 + PEE1(1,2)×1                                    | ECE(1,3)×1 [j→i] + PCE(1,3)×1 [j→i]                               | ECE(1,4)×1 [j→i] + PCP(1,4)×1 [j→i]                                |
| **2** | ENN(2,0)×1       | ENC(2,1)×1 + EEE(2,1)×2 + PEE0(2,1)×1 + PEE1(2,1)×1      | ENN(2,2)×1 + ENC(2,2)×2 [i→j] + ENC(2,2)×2 [j→i] + EEE(2,2)×4 + PEE0(2,2)×1 + PEE1(2,2)×1         | ECN(2,3)×1 + ECE(2,3)×2 [i→j] + PCE(2,3)×2 [i→j]                 | ECN(2,4)×1 + ECE(2,4)×2 [i→j] + PCP(2,4)×2 [i→j]                   |
| **3** | ECN(3,0)×1       | ECE(3,1)×1 + PCE(3,1)×1                                | ECN(3,2)×1 + ECE(3,2)×2 + PCE(3,2)×2                                         | ECC(3,3)×1                                             | ECC(3,4)×1                                              |
| **4** | ECN(4,0)×1       | ECE(4,1)×1 + PCP(4,1)×1                                | ECN(4,2)×1 + ECE(4,2)×2 + PCP(4,2)×2                                           | ECC(4,3)×1                                             | ECC(4,4)×1                                              |

### Function Calls Key

#### Basic Interaction Functions:
- `ElecNucNuc(qxq, rc, &ecoul, &fpair)`: Nucleus-nucleus Coulomb interaction
- `ElecNucElec(q, rc, eradius, &ecoul, &fpair, &force)`: Nucleus-electron interaction
- `ElecElecElec(rc, eradius1, eradius2, &ecoul, &fpair, &f1, &f2)`: Electron-electron interaction
- `ElecCoreElec(q, rc, eradius_core, eradius_elec, &ecoul, &fpair, &f1)`: Core-electron interaction
- `ElecCoreNuc(qxq, rc, eradius, &ecoul, &fpair)`: Core-nucleus interaction
- `ElecCoreCore(qxq, rc, eradius1, eradius2, &ecoul, &fpair)`: Core-core interaction

#### Pauli Repulsion Functions:
- `PauliElecElec(spin_same, rc, eradius1, eradius2, &epauli, &fpair, &f1, &f2)`: Pauli repulsion between electrons
- `PauliCoreElec(rc, eradius, &epauli, &fpair, &f1, A, B, C)`: Pauli repulsion for s-type pseudo-core (ECP_s)
- `PauliCorePElec(rc, eradius, &epauli, &fpair, &f1, A, B, C, D, E)`: Pauli repulsion for p-type pseudo-core (ECP_p)

## Cases

#### 1-0 Electron - Nucleus
```cpp
else if (abs(spin[i]) == 1 && spin[j] == 0) {
  e1rforce = 0.0;
  ElecNucElec(q[j], rc, eradius[i], &ecoul, &fpair, &e1rforce);
  e1rforce = spline * qqrd2e * e1rforce;
  erforce[i] += e1rforce;
  if (evflag && pressure_with_evirials_flag) {
    e1rvirial = eradius[i] * e1rforce;
    ev_tally_eff(i, i, nlocal, newton_pair, 0.0, e1rvirial);
  }
}
```

#### 1-1 Electron - Electron
```cpp
else if (abs(spin[i]) == 1 && abs(spin[j]) == 1) {
  e1rforce = e2rforce = 0.0;
  s_e1rforce = s_e2rforce = 0.0;

  ElecElecElec(rc, eradius[i], eradius[j], &ecoul, &fpair, &e1rforce, &e2rforce);
  PauliElecElec(spin[i] == spin[j], rc, eradius[i], eradius[j], &epauli, &s_fpair, &s_e1rforce, &s_e2rforce);

  epauli *= hhmss2e;
  s_fpair *= hhmss2e;

  e1rforce = spline * (qqrd2e * e1rforce + hhmss2e * s_e1rforce);
  erforce[i] += e1rforce;
  e2rforce = spline * (qqrd2e * e2rforce + hhmss2e * s_e2rforce);
  erforce[j] += e2rforce;

  if (evflag && pressure_with_evirials_flag) {
    e1rvirial = eradius[i] * e1rforce;
    e2rvirial = eradius[j] * e2rforce;
    ev_tally_eff(i, j, nlocal, newton_pair, 0.0, e1rvirial + e2rvirial);
  }
}
```

#### 1-2 Electron - Fixed-Core
```cpp
else if (abs(spin[i]) == 1 && spin[j] == 2) {
  e1rforce = e2rforce = 0.0;
  s_e1rforce = s_e2rforce = 0.0;

  ElecNucElec(q[j], rc, eradius[i], &ecoul, &fpair, &e2rforce);
  ElecElecElec(rc, eradius[j], eradius[i], &ecoul, &fpair, &e1rforce, &e2rforce);
  ElecElecElec(rc, eradius[j], eradius[i], &ecoul, &fpair, &e1rforce, &e2rforce);

  PauliElecElec(0, rc, eradius[j], eradius[i], &epauli, &s_fpair, &s_e1rforce, &s_e2rforce);
  PauliElecElec(1, rc, eradius[j], eradius[i], &epauli, &s_fpair, &s_e1rforce, &s_e2rforce);

  epauli *= hhmss2e;
  s_fpair *= hhmss2e;

  e2rforce = spline * (qqrd2e * e2rforce + hhmss2e * s_e2rforce);
  erforce[i] += e2rforce;
  if (evflag && pressure_with_evirials_flag) {
    e2rvirial = eradius[i] * e2rforce;
    ev_tally_eff(i, i, nlocal, newton_pair, 0.0, e2rvirial);
  }
}
```

#### 2-1 Fixed-Core - Electron
```cpp
else if (spin[i] == 2 && abs(spin[j]) == 1) {
  e1rforce = e2rforce = 0.0;
  s_e1rforce = s_e2rforce = 0.0;

  ElecNucElec(q[i], rc, eradius[j], &ecoul, &fpair, &e2rforce);
  ElecElecElec(rc, eradius[i], eradius[j], &ecoul, &fpair, &e1rforce, &e2rforce);
  ElecElecElec(rc, eradius[i], eradius[j], &ecoul, &fpair, &e1rforce, &e2rforce);

  PauliElecElec(0, rc, eradius[i], eradius[j], &epauli, &s_fpair, &s_e1rforce, &s_e2rforce);
  PauliElecElec(1, rc, eradius[i], eradius[j], &epauli, &s_fpair, &s_e1rforce, &s_e2rforce);

  epauli *= hhmss2e;
  s_fpair *= hhmss2e;

  e2rforce = spline * (qqrd2e * e2rforce + hhmss2e * s_e2rforce);
  erforce[j] += e2rforce;
  if (evflag && pressure_with_evirials_flag) {
    e2rvirial = eradius[j] * e2rforce;
    ev_tally_eff(j, j, nlocal, newton_pair, 0.0, e2rvirial);
  }
}
```

#### 2-2 Fixed-Core - Fixed-Core
```cpp
else if (spin[i] == 2 && spin[j] == 2) {
  e1rforce = e2rforce = 0.0;
  s_e1rforce = s_e2rforce = 0.0;
  double qxq = q[i] * q[j];

  ElecNucNuc(qxq, rc, &ecoul, &fpair);
  ElecNucElec(q[i], rc, eradius[j], &ecoul, &fpair, &e1rforce);
  ElecNucElec(q[i], rc, eradius[j], &ecoul, &fpair, &e1rforce);
  ElecNucElec(q[j], rc, eradius[i], &ecoul, &fpair, &e1rforce);
  ElecNucElec(q[j], rc, eradius[i], &ecoul, &fpair, &e1rforce);
  ElecElecElec(rc, eradius[i], eradius[j], &ecoul, &fpair, &e1rforce, &e2rforce);
  ElecElecElec(rc, eradius[i], eradius[j], &ecoul, &fpair, &e1rforce, &e2rforce);
  ElecElecElec(rc, eradius[i], eradius[j], &ecoul, &fpair, &e1rforce, &e2rforce);
  ElecElecElec(rc, eradius[i], eradius[j], &ecoul, &fpair, &e1rforce, &e2rforce);
  PauliElecElec(0, rc, eradius[i], eradius[j], &epauli, &s_fpair, &s_e1rforce, &s_e2rforce);
  PauliElecElec(1, rc, eradius[i], eradius[j], &epauli, &s_fpair, &s_e1rforce, &s_e2rforce);
  epauli *= 2;
  s_fpair *= 2;
  epauli *= hhmss2e;
  s_fpair *= hhmss2e;
}
```

#### 1-3 Electron - Pseudo-Core
```cpp
else if (abs(spin[i]) == 1 && spin[j] == 3) {
  e1rforce = ecp_e1rforce = 0.0;

  if (PAULI_CORE_D[ecp_type[jtype]] == 0.0 && PAULI_CORE_E[ecp_type[jtype]] == 0.0) {
    ElecCoreElec(q[j], rc, eradius[j], eradius[i], &ecoul, &fpair, &e1rforce);
    PauliCoreElec(rc, eradius[i], &ecp_epauli, &ecp_fpair, &ecp_e1rforce,
                  PAULI_CORE_A[ecp_type[jtype]],
                  PAULI_CORE_B[ecp_type[jtype]],
                  PAULI_CORE_C[ecp_type[jtype]]);
  } else {
    ElecCoreElec(q[j], rc, eradius[j], eradius[i], &ecoul, &fpair, &e1rforce);
    PauliCorePElec(rc, eradius[i], &ecp_epauli, &ecp_fpair, &ecp_e1rforce,
                   PAULI_CORE_A[ecp_type[jtype]],
                   PAULI_CORE_B[ecp_type[jtype]],
                   PAULI_CORE_C[ecp_type[jtype]],
                   PAULI_CORE_D[ecp_type[jtype]],
                   PAULI_CORE_E[ecp_type[jtype]]);
  }

  ecp_epauli *= h2e;
  ecp_fpair *= h2e;
  e1rforce = spline * (qqrd2e * e1rforce + h2e * ecp_e1rforce);
  erforce[i] += e1rforce;
  if (evflag && pressure_with_evirials_flag) {
    e1rvirial = eradius[i] * e1rforce;
    ev_tally_eff(i, i, nlocal, newton_pair, 0.0, e1rvirial);
  }
}
```

#### 2-3 Fixed-Core - Pseudo-Core
```cpp
else if (spin[i] == 2 && spin[j] == 3) {
  e2rforce = ecp_e2rforce = 0.0;

  if (PAULI_CORE_D[ecp_type[itype]] == 0.0 && PAULI_CORE_E[ecp_type[itype]] == 0.0) {
    double qxq = q[i] * q[j];
    ElecCoreNuc(qxq, rc, eradius[j], &ecoul, &fpair);
    ElecCoreElec(q[i], rc, eradius[i], eradius[j], &ecoul, &fpair, &e2rforce);
    ElecCoreElec(q[i], rc, eradius[i], eradius[j], &ecoul, &fpair, &e2rforce);
    PauliCoreElec(rc, eradius[j], &ecp_epauli, &ecp_fpair, &ecp_e2rforce,
                  PAULI_CORE_A[ecp_type[itype]],
                  PAULI_CORE_B[ecp_type[itype]],
                  PAULI_CORE_C[ecp_type[itype]]);
  } else {
    double qxq = q[i] * q[j];
    ElecCoreNuc(qxq, rc, eradius[j], &ecoul, &fpair);
    ElecCoreElec(q[i], rc, eradius[i], eradius[j], &ecoul, &fpair, &e2rforce);
    ElecCoreElec(q[i], rc, eradius[i], eradius[j], &ecoul, &fpair, &e2rforce);
    PauliCorePElec(rc, eradius[j], &ecp_epauli, &ecp_fpair, &ecp_e2rforce,
                   PAULI_CORE_A[ecp_type[itype]],
                   PAULI_CORE_B[ecp_type[itype]],
                   PAULI_CORE_C[ecp_type[itype]],
                   PAULI_CORE_D[ecp_type[itype]],
                   PAULI_CORE_E[ecp_type[itype]]);
  }

  ecp_epauli *= h2e;
  ecp_fpair *= h2e;
  e2rforce = spline * (qqrd2e * e2rforce + h2e * ecp_e2rforce);
  erforce[j] += e2rforce;
  if (evflag && pressure_with_evirials_flag) {
    e2rvirial = eradius[j] * e2rforce;
    ev_tally_eff(j, j, nlocal, newton_pair, 0.0, e2rvirial);
  }
}
```

#### 3-1 Pseudo-Core - Electron
```cpp
else if (spin[i] == 3 && abs(spin[j]) == 1) {
  e2rforce = ecp_e2rforce = 0.0;

  ElecCoreElec(q[i], rc, eradius[i], eradius[j], &ecoul, &fpair, &e2rforce);
  PauliCoreElec(rc, eradius[j], &ecp_epauli, &ecp_fpair, &ecp_e2rforce,
                PAULI_CORE_A[ecp_type[itype]],
                PAULI_CORE_B[ecp_type[itype]],
                PAULI_CORE_C[ecp_type[itype]]);

  ecp_epauli *= h2e;
  ecp_fpair *= h2e;
  e2rforce = spline * (qqrd2e * e2rforce + h2e * ecp_e2rforce);
  erforce[j] += e2rforce;
  if (evflag && pressure_with_evirials_flag) {
    e2rvirial = eradius[j] * e2rforce;
    ev_tally_eff(j, j, nlocal, newton_pair, 0.0, e2rvirial);
  }
}
```

#### 3-2 Pseudo-Core - Fixed-Core
```cpp
else if (spin[i] == 3 && spin[j] == 2) {
  e2rforce = ecp_e2rforce = 0.0;

  double qxq = q[i] * q[j];
  ElecCoreNuc(qxq, rc, eradius[j], &ecoul, &fpair);
  ElecCoreElec(q[i], rc, eradius[i], eradius[j], &ecoul, &fpair, &e2rforce);
  ElecCoreElec(q[i], rc, eradius[i], eradius[j], &ecoul, &fpair, &e2rforce);
  PauliCorePElec(rc, eradius[j], &ecp_epauli, &ecp_fpair, &ecp_e2rforce,
                 PAULI_CORE_A[ecp_type[itype]],
                 PAULI_CORE_B[ecp_type[itype]],
                 PAULI_CORE_C[ecp_type[itype]],
                 PAULI_CORE_D[ecp_type[itype]],
                 PAULI_CORE_E[ecp_type[itype]]);

  ecp_epauli *= h2e;
  ecp_fpair *= h2e;
  e2rforce = spline * (qqrd2e * e2rforce + h2e * ecp_e2rforce);
  erforce[j] += e2rforce;
  if (evflag && pressure_with_evirials_flag) {
    e2rvirial = eradius[j] * e2rforce;
    ev_tally_eff(j, j, nlocal, newton_pair, 0.0, e2rvirial);
  }
}
```
{{ ... }}
