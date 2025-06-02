#pragma once

/* -*- c++ -*- ----------------------------------------------------------
 LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
 https://www.lammps.org/, Sandia National Laboratories
 LAMMPS development team: developers@lammps.org

 Copyright (2003) Sandia Corporation.  Under the terms of Contract
 DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
 certain rights in this software.  This software is distributed under
 the GNU General Public License.

 See the README file in the top-level LAMMPS directory.
 ------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
 Contributing authors: Andres Jaramillo-Botero, Hai Xiao, Julius Su (Caltech)
------------------------------------------------------------------------- */

namespace LAMMPS_NS {

    #define PAULI_RE 0.9
    #define PAULI_RC 1.125
    #define PAULI_RHO -0.2
    
    #define ERF_TERMS1 12
    #define ERF_TERMS2 7
    #define DERF_TERMS 13
    
    // error arrays
    
    static constexpr double E1[] = {1.483110564084803581889448079057,
                                    -3.01071073386594942470731046311E-1,
                                    6.8994830689831566246603180718E-2,
                                    -1.3916271264722187682546525687E-2,
                                    2.420799522433463662891678239E-3,
                                    -3.65863968584808644649382577E-4,
                                    4.8620984432319048282887568E-5,
                                    -5.749256558035684835054215E-6,
                                    6.11324357843476469706758E-7,
                                    -5.8991015312958434390846E-8,
                                    5.207009092068648240455E-9,
                                    -4.23297587996554326810E-10,
                                    3.1881135066491749748E-11,
                                    -2.236155018832684273E-12,
                                    1.46732984799108492E-13,
                                    -9.044001985381747E-15,
                                    5.25481371547092E-16,
                                    -2.8874261222849E-17,
                                    1.504785187558E-18,
                                    -7.4572892821E-20,
                                    3.522563810E-21,
                                    -1.58944644E-22,
                                    6.864365E-24,
                                    -2.84257E-25,
                                    1.1306E-26,
                                    -4.33E-28,
                                    1.6E-29,
                                    -1.0E-30};
    
    static constexpr double E2[] = {1.077977852072383151168335910348,
                                    -2.6559890409148673372146500904E-2,
                                    -1.487073146698099509605046333E-3,
                                    -1.38040145414143859607708920E-4,
                                    -1.1280303332287491498507366E-5,
                                    -1.172869842743725224053739E-6,
                                    -1.03476150393304615537382E-7,
                                    -1.1899114085892438254447E-8,
                                    -1.016222544989498640476E-9,
                                    -1.37895716146965692169E-10,
                                    -9.369613033737303335E-12,
                                    -1.918809583959525349E-12,
                                    -3.7573017201993707E-14,
                                    -3.7053726026983357E-14,
                                    2.627565423490371E-15,
                                    -1.121322876437933E-15,
                                    1.84136028922538E-16,
                                    -4.9130256574886E-17,
                                    1.0704455167373E-17,
                                    -2.671893662405E-18,
                                    6.49326867976E-19,
                                    -1.65399353183E-19,
                                    4.2605626604E-20,
                                    -1.1255840765E-20,
                                    3.025617448E-21,
                                    -8.29042146E-22,
                                    2.31049558E-22,
                                    -6.5469511E-23,
                                    1.8842314E-23,
                                    -5.504341E-24,
                                    1.630950E-24,
                                    -4.89860E-25,
                                    1.49054E-25,
                                    -4.5922E-26,
                                    1.4318E-26,
                                    -4.516E-27,
                                    1.440E-27,
                                    -4.64E-28,
                                    1.51E-28,
                                    -5.0E-29,
                                    1.7E-29,
                                    -6.0E-30,
                                    2.0E-30,
                                    -1.0E-30};
    
    static constexpr double DE1[] = {-0.689379974848418501361491576718,
                                     0.295939056851161774752959335568,
                                     -0.087237828075228616420029484096,
                                     0.019959734091835509766546612696,
                                     -0.003740200486895490324750329974,
                                     0.000593337912367800463413186784,
                                     -0.000081560801047403878256504204,
                                     9.886099179971884018535968E-6,
                                     -1.071209234904290565745194E-6,
                                     1.0490945447626050322784E-7,
                                     -9.370959271038746709966E-9,
                                     7.6927263488753841874E-10,
                                     -5.8412335114551520146E-11,
                                     4.125393291736424788E-12,
                                     -2.72304624901729048E-13,
                                     1.6869717361387012E-14,
                                     -9.84565340276638E-16,
                                     5.4313471880068E-17,
                                     -2.840458699772E-18,
                                     1.4120512798E-19,
                                     -6.688772574E-21,
                                     3.0257558E-22,
                                     -1.3097526E-23,
                                     5.4352E-25,
                                     -2.1704E-26,
                                     8.32E-28,
                                     -5.4E-29};
    
    static constexpr double DE2[] = {0.717710208167480928473053690384,
                                     -0.379868973985143305103199928808,
                                     0.125832094465157378967135019248,
                                     -0.030917661684228839423081992424,
                                     0.006073689914144320367855343072,
                                     -0.000996057789064916825079352632,
                                     0.000140310790466315733723475232,
                                     -0.000017328176496070286001302184,
                                     1.90540194670935746397168e-6,
                                     -1.8882873760163694937908e-7,
                                     1.703176613666840587056e-8,
                                     -1.40955218086201517976e-9,
                                     1.0776816914256065828e-10,
                                     -7.656138112778696256e-12,
                                     5.07943557413613792e-13,
                                     -3.1608615530282912e-14,
                                     1.852036572003432e-15,
                                     -1.02524641430496e-16,
                                     5.37852808112e-18,
                                     -2.68128238704e-19,
                                     1.273321788e-20,
                                     -5.77335744e-22,
                                     2.504352e-23,
                                     -1.0446e-24,
                                     4.16e-26,
                                     -2.808e-27};
    
    // inline functions for performance
    
    // ----------------------------------------------------------------------
    
    static inline double ipoly02(double x){
      // P(x) in the range x > 2
      int i;
      double b0, b1, b2;
      b1 = 0.0;
      b0 = 0.0;
      x *= 2;
      for (i = ERF_TERMS2; i >= 0; i--) {
        b2 = b1;
        b1 = b0;
        b0 = x * b1 - b2 + E2[i];
      }
      return 0.5 * (b0 - b2);
    }
    
    // ----------------------------------------------------------------------
    
    static inline double ipoly1(double x){
      // First derivative P'(x) in the range x < 2
      int i;
      double b0, b1, b2;
    
      b1 = 0.0;
      b0 = 0.0;
      x *= 2;
      for (i = DERF_TERMS; i >= 0; i--) {
        b2 = b1;
        b1 = b0;
        b0 = x * b1 - b2 + DE1[i];
      }
      return 0.5 * (b0 - b2);
    }
    
    // ----------------------------------------------------------------------
    
    static inline double ipoly01(double x){
      // P(x) in the range x < 2
    
      int i;
      double b0, b1, b2;
      b1 = 0.0;
      b0 = 0.0;
      x *= 2;
      for (i = ERF_TERMS1; i >= 0; i--) {
        b2 = b1;
        b1 = b0;
        b0 = x * b1 - b2 + E1[i];
      }
      return 0.5 * (b0 - b2);
    }
    
    // ----------------------------------------------------------------------
    
    static inline double ierfoverx1(double x, double *df){
      // Computes Erf(x)/x and its first derivative
    
      double t, f;
      double x2;    // x squared
      double exp_term, recip_x;
    
      if (x < 2.0) {
        // erf(x) = x * y(t)     
        // t = 2 * (x/2)^2 - 1.  
        t = 0.5 * x * x - 1;
        f = ipoly01(t);
        *df = ipoly1(t) * x;
      } else {
        // erf(x) = 1 - exp(-x^2)/x * y(t)
        // t = (10.5 - x^2) / (2.5 + x^2)
        x2 = x * x;
        t = (10.5 - x2) / (2.5 + x2);
        exp_term = exp(-x2);
        recip_x = 1.0 / x;
        f = 1.0 / x - (exp_term / x2) * ipoly02(t);
        *df = (1.12837916709551257389615890312 * exp_term - f) * recip_x;
      }
      return f;
    }
    
    // ----------------------------------------------------------------------
    
    static inline void KinElec(double radius, double *eke, double *frc){
      *eke += 1.5 / (radius * radius);
      *frc += 3.0 / (radius * radius * radius);
    }
    
    // ----------------------------------------------------------------------
    
    static inline void ElecNucNuc(double q, double rc, double *ecoul, double *frc){
      *ecoul += q / rc;
      *frc += q / (rc * rc);
    }
    
    /* ---------------------------------------------------------------------- */
    
    static inline void ElecNucElec(double q, double rc, double re1, double *ecoul, double *frc, double *fre1){
      double a, arc;
      double coeff_a;
    
      coeff_a = 1.4142135623730951; // sqrt(2)
    
      // E = -Z/r Erf(a r / re) 
      a = coeff_a / re1;
      arc = a * rc;
    
      // Interaction between nuclear point charge and Gaussian electron
      double E, dEdr, dEdr1, f, df;
    
      f = ierfoverx1(arc, &df);
      dEdr = q * a * a * df;
      dEdr1 = -q * (a / re1) * (f + arc * df);
      E = -q * a * f;
    
      *ecoul += E;
      *frc += dEdr;
      *fre1 += dEdr1;
    }
    
    // ----------------------------------------------------------------------
    
    static inline void ElecElecElec(double rc, double re1, double re2, double *ecoul, double *frc,double *fre1, double *fre2){
      double a, arc, re, fre;
      double coeff_a;
    
      coeff_a = 1.4142135623730951;  // sqrt(2)
      re = sqrt(re1 * re1 + re2 * re2);
      a = coeff_a / re;  
      arc = a * rc;       
    
    //   V_elecelec = E * F                             
    //   where E = -Z/r Erf(a r / re)                   
    //          F = (1 - (b s + c s^2) exp(-d s^2))      
    //   and s = r / re                                 
    
      double E, dEdr, dEdr1, dEdr2, f, df;
    
      f = ierfoverx1(arc, &df);
      dEdr = -a * a * df;
      fre = a * (f + arc * df) / (re * re);
      dEdr1 = fre * re1;
      dEdr2 = fre * re2;
    
      E = a * f;
    
      *ecoul += E;
      *frc += dEdr;
      *fre1 += dEdr1;
      *fre2 += dEdr2;
    }
    
    // ----------------------------------------------------------------------
    
    static inline void ElecCoreNuc(double q, double rc, double re1, double *ecoul, double *frc){
      double a, arc;
      double coeff_a;
      double E, dEdr, df, f;
    
      coeff_a = 1.4142135623730951; /* sqrt(2) */
      a = coeff_a / re1;
      arc = a * rc;
    
      f = ierfoverx1(arc, &df);
      dEdr = -q * a * a * df;
      E = q * a * f;
    
      *ecoul += E;
      *frc += dEdr;
    }
    
    // ----------------------------------------------------------------------
    
    static inline void ElecCoreCore(double q, double rc, double re1, double re2, double *ecoul, double *frc){
      double a, arc, re;
      double coeff_a;
      double E, dEdr, f, df;
    
      coeff_a = 1.4142135623730951;
    
      re = sqrt(re1 * re1 + re2 * re2);
      a = coeff_a / re;
      arc = a * rc;
    
      f = ierfoverx1(arc, &df);
      dEdr = -q * a * a * df;
      E = q * a * f;
    
      *ecoul += E;
      *frc += dEdr;
    }
    
    // ----------------------------------------------------------------------
    
    static inline void ElecCoreElec(double q, double rc, double re1, double re2, double *ecoul, double *frc, double *fre2){
      double a, arc, re;
      double coeff_a;
      double E, dEdr, dEdr2, f, df, fre;
    
      coeff_a = 1.4142135623730951;
    
      
    //   re1: core size
    //   re2: electron size
    //   re3: size of the core, obtained from the electron density function rho(r) of core
    //   e.g. rho(r) = a1*exp(-((r)/b1)^2), a1 =157.9, b1 = 0.1441 -> re3 = 0.1441 for Si4+
    
      re = sqrt(re1 * re1 + re2 * re2);
    
      a = coeff_a / re;
      arc = a * rc;
    
      f = ierfoverx1(arc, &df);
      E = -q * a * f;
      dEdr = -q * a * df * a;
      fre = q * a * (f + arc * df) / (re * re);
      dEdr2 = fre * re2;
    
      *ecoul += E;
      *frc -= dEdr;
      *fre2 -= dEdr2;
    }
    
    
    /*
     * @brief Calculates the Pauli repulsion energy and forces between two electrons.
     *
     * This function implements the electron-electron Pauli interaction term as used
     * in the eFF (electron Force Field). The interaction depends on the distance
     * between electrons, their sizes, and their relative spin orientation.
     *
     * The input parameters `re1` (size of electron 1), `re2` (size of electron 2),
     * and `rc` (distance between electrons) are first scaled using global macros:
     * \f[ re_1' = re1 \cdot \text{PAULI\_RE} \f]
     * \f[ re_2' = re2 \cdot \text{PAULI\_RE} \f]
     * \f[ r_c' = rc \cdot \text{PAULI\_RC} \f]
     *
     * The Pauli repulsion energy \f$E_P\f$ is then computed as \f$E_P = T_{eff} \cdot O(S)\f$, where:
     * 1. \f$S\f$ is the overlap integral between the scaled Gaussian electron wavefunctions:
     *    \f[ S(r_c', re_1', re_2') = \left( \frac{2 re_1' re_2'}{re_1'^2 + re_2'^2} \right)^{3/2} \exp\left(-\frac{r_c'^2}{re_1'^2 + re_2'^2}\right) \f]
     * 2. \f$T_{eff}\f$ is an effective kinetic energy term:
     *    \f[ T_{eff}(r_c', re_1', re_2') = \frac{3}{2}\left(\frac{1}{re_1'^2} + \frac{1}{re_2'^2}\right) - \frac{2(3(re_1'^2+re_2'^2) - 2r_c'^2)}{(re_1'^2+re_2'^2)^2} \f]
     * 3. \f$O(S)\f$ is a spin-dependent term, where \f$\rho = \text{PAULI\_RHO}\f$:
     *    - For same-spin electrons (`samespin == 1`):
     *      \f[ O(S) = \frac{S^2}{1 - S^2} + (1 - \rho) \frac{S^2}{1 + S^2} \f]
     *    - For opposite-spin electrons (`samespin == 0`):
     *      \f[ O(S) = -\rho \frac{S^2}{1 + S^2} \f]
     *
     * The calculated energy \f$E_P\f$ is added to `*epauli`.
     * The function also computes the negative derivatives of \f$E_P\f$ with respect to the
     * original (unscaled) \f$rc, re1, re2\f$ and adds them to `*frc, *fre1, *fre2` respectively.
     *
     * @param samespin Integer flag: 1 if electrons have the same spin, 0 otherwise.
     * @param rc Distance between the two electrons (unscaled).
     * @param re1 Size of the first electron (unscaled).
     * @param re2 Size of the second electron (unscaled).
     * @param epauli Pointer to which the calculated Pauli energy \f$E_P\f$ is added.
     * @param frc Pointer to which the force component \f$-\frac{\partial E_P}{\partial rc}\f$ is added.
     * @param fre1 Pointer to which the force component \f$-\frac{\partial E_P}{\partial re1}\f$ (force on electron 1 size) is added.
     * @param fre2 Pointer to which the force component \f$-\frac{\partial E_P}{\partial re2}\f$ (force on electron 2 size) is added.
     */
    static inline void PauliElecElec(int samespin, double rc, double re1, double re2, double *epauli, double *frc, double *fre1, double *fre2){
      double ree, rem;
      double S, t1, t2, tt;
      double dSdr1, dSdr2, dSdr;
      double dTdr1, dTdr2, dTdr;
      double O, dOdS, ratio;
    
      re1 *= PAULI_RE;
      re2 *= PAULI_RE;
      rc *= PAULI_RC;
      ree = re1 * re1 + re2 * re2;
      rem = re1 * re1 - re2 * re2;
    
      S = (2.82842712474619 / pow((re2 / re1 + re1 / re2), 1.5)) * exp(-rc * rc / ree);
    
      t1 = 1.5 * (1 / (re1 * re1) + 1 / (re2 * re2));
      t2 = 2.0 * (3 * ree - 2 * rc * rc) / (ree * ree);
      tt = t1 - t2;
    
      dSdr1 = (-1.5 / re1) * (rem / ree) + 2 * re1 * rc * rc / (ree * ree);
      dSdr2 = (1.5 / re2) * (rem / ree) + 2 * re2 * rc * rc / (ree * ree);
      dSdr  = -2 * rc / ree;
      dTdr1 = -3 / (re1 * re1 * re1) - 12 * re1 / (ree * ree) +    8 * re1 * (-2 * rc * rc + 3 * ree) / (ree * ree * ree);
      dTdr2 = -3 / (re2 * re2 * re2) - 12 * re2 / (ree * ree) +    8 * re2 * (-2 * rc * rc + 3 * ree) / (ree * ree * ree);
      dTdr  = 8 * rc / (ree * ree);
    
      if (samespin == 1) {
        O    = S * S /  (1.0 - S * S) + (1 - PAULI_RHO) * S * S / (1.0 + S * S);
        dOdS = 2 * S / ((1.0 - S * S) * (1.0 - S * S) ) + (1 - PAULI_RHO) * 2 * S / ((1.0 + S * S) * (1.0 + S * S));
      } else {
        O    = -PAULI_RHO * S * S / (1.0 + S * S);
        dOdS = -PAULI_RHO * 2 * S / ((1.0 + S * S) * (1.0 + S * S));
      }
    
      ratio = tt * dOdS * S;
      *fre1 -= PAULI_RE * (dTdr1 * O + ratio * dSdr1);
      *fre2 -= PAULI_RE * (dTdr2 * O + ratio * dSdr2);
      *frc -= PAULI_RC * (dTdr * O + ratio * dSdr);
      *epauli += tt * O;
    }
    
    
    /*
     * @brief Calculates the ECP Pauli repulsion energy and its derivatives for interactions  between a pseudo-core and an s-type valence electron.
     *
     * The implemented formula for energy (E) is:
     * \f[
     * E = \text{PAULI\_CORE\_A} \cdot \exp\left(-\frac{\text{PAULI\_CORE\_B} \cdot r_c^2}{s^2 + \text{PAULI\_CORE\_C}}\right)
     * \f]
     * The function calculates the energy $E$, its derivative with respect to $r_c$ ($\frac{\partial E}{\partial r_c}$),
     * and its derivative with respect to $s$ ($\frac{\partial E}{\partial s}$).
     *
     * @param rc Distance between the pseudo-core and the valence electron ($r_c$).
     * @param re2 Size of the valence electron ($s$).
     * @param epauli Pointer to which the calculated Pauli energy $E$ is added.
     * @param frc Pointer to which $-\frac{\partial E}{\partial r_c}$ is added.
     * @param fre2 Pointer to which $-\frac{\partial E}{\partial s}$ is added.
     * @param PAULI_CORE_A Parameter $P_A$ (pseudo-core wave function amplitude).
     * @param PAULI_CORE_B Parameter $P_B$ (pseudo-core wavefunction decay factor).
     * @param PAULI_CORE_C Parameter $P_C$ (square of effective pseudo-core particle size).
     */
    static inline void PauliCoreElec(double rc, double re2, double *epauli, double *frc, double *fre2, double PAULI_CORE_A, double PAULI_CORE_B, double PAULI_CORE_C) {
      double E, dEdrc, dEdre2, rcsq, ssq;
      rcsq = rc * rc;
      ssq = re2 * re2;
      E = PAULI_CORE_A * exp((-PAULI_CORE_B * rcsq) / (ssq + PAULI_CORE_C));
      dEdrc  = -2 * PAULI_CORE_A * PAULI_CORE_B * rc *          exp(-PAULI_CORE_B * rcsq / (ssq + PAULI_CORE_C)) /  (ssq + PAULI_CORE_C);
      dEdre2 =  2 * PAULI_CORE_A * PAULI_CORE_B * re2 * rcsq *  exp(-PAULI_CORE_B * rcsq / (ssq + PAULI_CORE_C)) /  ((PAULI_CORE_C + ssq) * (PAULI_CORE_C + ssq));
      *epauli += E;
      *frc -= dEdrc;
      *fre2 -= dEdre2;
    }
    
    /*
     * @brief Calculates the ECP Pauli repulsion energy and its derivatives for interactions between a pseudo-core and a p-type valence electron (represented by an s-type Gaussian).
     *
     * The implemented formula for energy (E) is:
     * \f[
     * E = \text{PAULI\_CORE\_P\_A} \left( \frac{2}{\frac{\text{PAULI\_CORE\_P\_B}}{s} + \frac{s}{\text{PAULI\_CORE\_P\_B}}} \right)^5 (r_c - \text{PAULI\_CORE\_P\_C} s)^2 \exp\left( -\frac{\text{PAULI\_CORE\_P\_D}(r_c - \text{PAULI\_CORE\_P\_C} s)^2}{\text{PAULI\_CORE\_P\_E} + s^2} \right)
     * \f]
     * Note: The prefactor $\left( \frac{2}{\frac{\text{PAULI\_CORE\_P\_B}}{s} + \frac{s}{\text{PAULI\_CORE\_P\_B}}} \right)$ is raised to the power of 5 in this implementation,
     * which differs from some documentation (e.g., ECP.md) where it might be to the power of 1.
     * The function calculates the energy $E$, its derivative with respect to $r_c$ ($\frac{\partial E}{\partial r_c}$),
     * and its derivative with respect to $s$ ($\frac{\partial E}{\partial s}$).
     *
     * @param rc Distance between the pseudo-core and the valence electron ($r_c$).
     * @param re2 Size of the valence electron ($s$).
     * @param epauli Pointer to which the calculated Pauli energy $E$ is added.
     * @param frc Pointer to which $-\frac{\partial E}{\partial r_c}$ is added.
     * @param fre2 Pointer to which $-\frac{\partial E}{\partial s}$ is added.
     * @param PAULI_CORE_P_A Parameter $P_A$ (pseudo-core wave function amplitude).
     * @param PAULI_CORE_P_B Parameter $P_B$ (second effective size for overlap amplitude).
     * @param PAULI_CORE_P_C Parameter $P_C$ (off-center measure).
     * @param PAULI_CORE_P_D Parameter $P_D$ (pseudo-core wavefunction decay factor).
     * @param PAULI_CORE_P_E Parameter $P_E$ (square of effective pseudo-core particle size).
     */
    static inline void PauliCorePElec(double rc, double re2, double *epauli, double *frc, double *fre2, double PAULI_CORE_P_A, double PAULI_CORE_P_B, double PAULI_CORE_P_C, double PAULI_CORE_P_D, double PAULI_CORE_P_E ){
      double E, dEdrc, dEdre2;
    
      E = PAULI_CORE_P_A * pow((2.0 / (PAULI_CORE_P_B / re2 + re2 / PAULI_CORE_P_B)), 5.0) * pow((rc - PAULI_CORE_P_C * re2), 2.0) *  exp(-PAULI_CORE_P_D * pow((rc - PAULI_CORE_P_C * re2), 2.0) / (PAULI_CORE_P_E + re2 * re2));
    
      dEdrc = PAULI_CORE_P_A * pow((2.0 / (PAULI_CORE_P_B / re2 + re2 / PAULI_CORE_P_B)), 5.0) * 2.0 * (rc - PAULI_CORE_P_C * re2) * exp(-PAULI_CORE_P_D * pow((rc - PAULI_CORE_P_C * re2), 2.0) /   (PAULI_CORE_P_E + re2 * re2)) 
            + E * (-PAULI_CORE_P_D * 2.0 * (rc - PAULI_CORE_P_C * re2) / (PAULI_CORE_P_E + re2 * re2));
    
      dEdre2 = E *
              (-5.0 / (PAULI_CORE_P_B / re2 + re2 / PAULI_CORE_P_B) *
               (-PAULI_CORE_P_B / (re2 * re2) + 1.0 / PAULI_CORE_P_B)) +
          PAULI_CORE_P_A * pow((2.0 / (PAULI_CORE_P_B / re2 + re2 / PAULI_CORE_P_B)), 5.0) * 2.0 *
              (rc - PAULI_CORE_P_C * re2) * (-PAULI_CORE_P_C) *
              exp(-PAULI_CORE_P_D * pow((rc - PAULI_CORE_P_C * re2), 2.0) /
                  (PAULI_CORE_P_E + re2 * re2)) 
          +
          E *
              (2.0 * PAULI_CORE_P_D * (rc - PAULI_CORE_P_C * re2) *
               (PAULI_CORE_P_C * PAULI_CORE_P_E + rc * re2) / pow((PAULI_CORE_P_E + re2 * re2), 2.0));
    
      *epauli += E;
      *frc -= dEdrc;
      *fre2 -= dEdre2;
    }
    
    /**
     * @brief Calculates the ECP Pauli repulsion energy and its derivatives for interactions
     *        between a pseudo-core and an s-type valence electron. (Optimized)
     *
     * The implemented formula for energy (E_calc) is:
     * \f[
     * E_{calc} = \text{paramA} \cdot \exp\left(-\frac{\text{paramB} \cdot rc^2}{re2^2 + \text{paramC}}\right)
     * \f]
     * The function calculates the energy $E_{calc}$, its derivative with respect to $rc$ ($\frac{\partial E_{calc}}{\partial rc}$),
     * and its derivative with respect to $re2$ ($\frac{\partial E_{calc}}{\partial re2}$).
     *
     * @param rc Distance between the pseudo-core and the valence electron ($rc$).
     * @param re2 Size of the valence electron ($s$ in the formula, passed as re2).
     * @param epauli Pointer to which the calculated Pauli energy $E_{calc}$ is added.
     * @param frc Pointer to which $-\frac{\partial E_{calc}}{\partial rc}$ is added.
     * @param fre2 Pointer to which $-\frac{\partial E_{calc}}{\partial re2}$ is added.
     * @param paramA Parameter $P_A$ (pseudo-core wave function amplitude).
     * @param paramB Parameter $P_B$ (pseudo-core wavefunction decay factor).
     * @param paramC Parameter $P_C$ (square of effective pseudo-core particle size).
     */
    static inline void PauliCoreElec_optGemmini(double rc, double re2, double *epauli, double *frc, double *fre2,
        double paramA, double paramB, double paramC) {
    // Pre-calculate squares and common terms
    double rc_sq = rc * rc;
    double re2_sq = re2 * re2;

    double den = re2_sq + paramC;
    // It's generally assumed den > 0 (paramC is square of a size, re2_sq is square of a size)
    // Add a check or ensure params guarantee this if necessary.
    // if (den == 0.0) { /* handle error or return */ return; } 
    double inv_den = 1.0 / den;

    double exp_arg = -paramB * rc_sq * inv_den;
    double exp_val = exp(exp_arg);

    double E_calc = paramA * exp_val;

    // Derivatives can reuse E_calc and other pre-calculated terms
    // dE/d(rc) = A * exp(...) * (-B * 2*rc / den)
    //          = E_calc * (-2 * paramB * rc * inv_den)
    double dEdrc = E_calc * (-2.0 * paramB * rc * inv_den);

    // dE/d(re2) = A * exp(...) * (-B * rc_sq * (-1) * (2*re2) / den^2)
    //           = E_calc * (2 * paramB * re2 * rc_sq * inv_den * inv_den)
    double dEdre2 = E_calc * (2.0 * paramB * re2 * rc_sq * inv_den * inv_den);

    *epauli += E_calc;
    *frc -= dEdrc;    // Add -dE/drc
    *fre2 -= dEdre2;  // Add -dE/dre2
    }

   
    static inline void PauliCoreElec_optDeepSeek(double rc, double re2, double *epauli, double *frc, double *fre2, double A, double B, double C) {
        double rcsq  = rc * rc;
        double ssq   = re2 * re2;
        double denom = ssq + C;
        double inv_denom = 1.0 / denom;
        double arg = -B * rcsq * inv_denom;
        double exp_val = exp(arg);
        double coef = A * B * exp_val;
        
        double E = A * exp_val;
        double dEdrc = -2.0 * coef * rc * inv_denom;
        double dEdre2 = 2.0 * coef * re2 * rcsq * inv_denom * inv_denom;
        
        *epauli += E;
        *frc -= dEdrc;
        *fre2 -= dEdre2;
    }
    
    static inline void PauliCorePElec_optDeepSeek(double rc, double re2, double *epauli, double *frc, double *fre2, double A, double B, double C, double D, double E_param) {
        double re2_2       = re2 * re2;
        double B2          = B * B;
        double denom_ratio = B2 + re2_2;
        double ratio       = 2.0 * re2 * B / denom_ratio;
        
        // Compute ratio^5 efficiently
        double ratio2 = ratio  * ratio;
        double ratio4 = ratio2 * ratio2;
        double ratio5 = ratio4 * ratio;
        
        double delta    = rc - C * re2;
        double delta2   = delta * delta;
        double denom2   = E_param + re2_2;
        double inv_denom2 = 1.0 / denom2;
        double exp_arg  = -D * delta2 * inv_denom2;
        double exp_val  = exp(exp_arg);
        
        double term1_rc = A * ratio5 * exp_val;  // A * ratio^5 * exp_val
        double E_val = term1_rc * delta2;        // Full energy value
        
        // Precompute common terms for derivatives
        double U = D * delta2 * inv_denom2;
        double dEdrc_val = 2.0 * delta * term1_rc * (1.0 - U);
        
        // dEdre2 components
        double term1 = (denom_ratio != 0.0 && re2 != 0.0) ?   E_val * 5.0 * (B2 - re2_2) / (denom_ratio * re2) : 0.0;
        double term2 = term1_rc * 2.0 * (-C) * delta;
        double term3 = E_val    * 2.0 *   D  * delta * (C * denom2 + re2 * delta) * (inv_denom2 * inv_denom2);
        double dEdre2_val = term1 + term2 + term3;
        
        *epauli += E_val;
        *frc -= dEdrc_val;
        *fre2 -= dEdre2_val;
    }

    /**
     * @brief Calculates the ECP Pauli repulsion energy and its derivatives for interactions
     *        between a pseudo-core and a p-type valence electron (represented by an s-type Gaussian). (Optimized)
     *
     * The implemented formula for energy (E_calc) is:
     * \f[
     * E_{calc} = \text{paramA} \left( \frac{2}{\frac{\text{paramB}}{re2} + \frac{re2}{\text{paramB}}} \right)^5 (rc - \text{paramC} re2)^2 \exp\left( -\frac{\text{paramD}(rc - \text{paramC} re2)^2}{\text{paramE} + re2^2} \right)
     * \f]
     * The function calculates the energy $E_{calc}$, its derivative with respect to $rc$ ($\frac{\partial E_{calc}}{\partial rc}$),
     * and its derivative with respect to $re2$ ($\frac{\partial E_{calc}}{\partial re2}$).
     *
     * @param rc Distance between the pseudo-core and the valence electron ($rc$).
     * @param re2 Size of the valence electron ($s$ in the formula, passed as re2).
     * @param epauli Pointer to which the calculated Pauli energy $E_{calc}$ is added.
     * @param frc Pointer to which $-\frac{\partial E_{calc}}{\partial rc}$ is added.
     * @param fre2 Pointer to which $-\frac{\partial E_{calc}}{\partial re2}$ is added.
     * @param paramA Parameter $P_A$ (pseudo-core wave function amplitude).
     * @param paramB Parameter $P_B$ (second effective size for overlap amplitude).
     * @param paramC Parameter $P_C$ (off-center measure).
     * @param paramD Parameter $P_D$ (pseudo-core wavefunction decay factor).
     * @param paramE Parameter $P_E$ (square of effective pseudo-core particle size).
     */
    static inline void PauliCorePElec_optGemini(double rc, double re2, double *epauli, double *frc, double *fre2,
        double paramA, double paramB, double paramC, double paramD, double paramE_const) {
// It's assumed re2 > 0 and paramB > 0. Add checks if necessary.
// if (re2 <= 0.0 || paramB <= 0.0) { /* handle error or return */ return; }

double re2_sq = re2 * re2;
double inv_re2 = 1.0 / re2;

// Prefactor term: (2 / (paramB/re2 + re2/paramB))^5
double B_div_s = paramB * inv_re2;
double s_div_B = re2 / paramB; // or re2 * inv_paramB if paramB is constant over many calls
double sum_ratios = B_div_s + s_div_B;
// if (sum_ratios == 0.0) { /* handle error or return */ return; }
double inv_sum_ratios = 1.0 / sum_ratios;
double prefactor_base = 2.0 * inv_sum_ratios;

double pb_sq = prefactor_base * prefactor_base;
double pb_p4 = pb_sq * pb_sq;
double prefactor_pow5 = pb_p4 * prefactor_base;

// Effective distance term: (rc - paramC * re2)
double R_eff = rc - paramC * re2;
double R_eff_sq = R_eff * R_eff;

// Exponential term
double exp_den = paramE_const + re2_sq;
// if (exp_den == 0.0) { /* handle error or return */ return; }
double inv_exp_den = 1.0 / exp_den;
double exp_arg = -paramD * R_eff_sq * inv_exp_den;
double exp_val = exp(exp_arg);

// Energy
double E_calc = paramA * prefactor_pow5 * R_eff_sq * exp_val;

// Derivatives
// Common factor for some derivative terms
double term_A_P5_Exp = paramA * prefactor_pow5 * exp_val;

// dE/drc = A*P^5 * exp(...) * d(R_eff_sq)/drc + A*P^5*R_eff_sq * d(exp(...))/drc
// d(R_eff_sq)/drc = 2*R_eff
// d(exp(...))/drc = exp(...) * (-D * d(R_eff_sq)/drc / exp_den)
//                 = exp(...) * (-D * 2*R_eff / exp_den)
// dEdrc = term_A_P5_Exp * 2.0 * R_eff + (paramA * prefactor_pow5 * R_eff_sq) * exp_val * (-paramD * 2.0 * R_eff * inv_exp_den)
//       = term_A_P5_Exp * 2.0 * R_eff + E_calc * (-paramD * 2.0 * R_eff * inv_exp_den)
// This form is robust if R_eff = 0 (making E_calc = 0 and dEdrc = 0)
double dEdrc = term_A_P5_Exp * 2.0 * R_eff + E_calc * (-2.0 * paramD * R_eff * inv_exp_den);

// dE/dre2 has three parts from differentiating Prefactor, R_eff_sq, and Exp_val w.r.t re2
// Part 1: d(Prefactor_pow5)/dre2
// d(sum_ratios)/dre2 = -paramB/(re2*re2) + 1.0/paramB
double d_sum_ratios_ds = -paramB * inv_re2 * inv_re2 + (1.0 / paramB);
// d(Prefactor_base)/dre2 = d(2.0/sum_ratios)/dre2 = -2.0 * inv_sum_ratios^2 * d_sum_ratios_ds
// d(Prefactor_pow5)/dre2 = 5 * Prefactor_base^4 * d(Prefactor_base)/dre2
//                        = 5 * (Prefactor_pow5 / Prefactor_base) * (-2.0 * inv_sum_ratios^2 * d_sum_ratios_ds)
//                        = Prefactor_pow5 * 5 * (sum_ratios/2.0) * (-2.0 * inv_sum_ratios^2 * d_sum_ratios_ds)
//                        = Prefactor_pow5 * (-5.0 * inv_sum_ratios * d_sum_ratios_ds)
double dEdre2_term1 = E_calc * (-5.0 * inv_sum_ratios * d_sum_ratios_ds);

// Part 2: d(R_eff_sq)/dre2
// d(R_eff_sq)/dre2 = 2 * R_eff * (-paramC)
double dEdre2_term2 = term_A_P5_Exp * 2.0 * R_eff * (-paramC);

// Part 3: d(Exp_val)/dre2
// d(Exp_arg)/dre2 = -paramD * [ d(R_eff_sq)/dre2 * inv_exp_den + R_eff_sq * d(inv_exp_den)/dre2 ]
// d(inv_exp_den)/dre2 = -1.0 * inv_exp_den^2 * (2.0*re2)
// d(Exp_arg)/dre2 = -paramD * [ (2*R_eff*(-paramC))*inv_exp_den - R_eff_sq*2*re2*inv_exp_den^2 ]
//                 = -paramD * 2 * R_eff * inv_exp_den^2 * [ (-paramC*exp_den) - (R_eff*re2) ]
//                 =  paramD * 2 * R_eff * inv_exp_den^2 * [ paramC*exp_den + R_eff*re2 ]
//                 =  paramD * 2 * R_eff * inv_exp_den^2 * [ paramC*(paramE_const+re2_sq) + (rc-paramC*re2)*re2 ]
//                 =  paramD * 2 * R_eff * inv_exp_den^2 * [ paramC*paramE_const + paramC*re2_sq + rc*re2 - paramC*re2_sq ]
//                 =  paramD * 2 * R_eff * inv_exp_den^2 * [ paramC*paramE_const + rc*re2 ]
// d(Exp_val)/dre2 = Exp_val * d(Exp_arg)/dre2
// So this term is E_calc * d(Exp_arg)/dre2 (if A_p*Prefactor_pow5*R_eff_sq is not zero)
// Original code: E_calc * (2.0 * D * R_eff * (C*E_const + rc*re2) / (exp_den^2))
// This matches my derivation for d(Exp_arg)/dre2.
double dEdre2_term3_factor_num = 2.0 * paramD * R_eff * (paramC * paramE_const + rc * re2);
double dEdre2_term3 = E_calc * (dEdre2_term3_factor_num * inv_exp_den * inv_exp_den);

double dEdre2 = dEdre2_term1 + dEdre2_term2 + dEdre2_term3;

*epauli += E_calc;
*frc -= dEdrc;
*fre2 -= dEdre2;
}




    // ----------------------------------------------------------------------
    
    static inline void RForce(double dx, double dy, double dz, double rc, double force, double *fx, double *fy, double *fz){
      force /= rc;
      *fx = force * dx;
      *fy = force * dy;
      *fz = force * dz;
    }
    
    // ----------------------------------------------------------------------
    
    static inline void SmallRForce(double dx, double dy, double dz, double rc, double force, double *fx, double *fy, double *fz){
      /* Handles case where rc is small to avoid division by zero */
      if (rc > 0.000001) {
        force /= rc;
        *fx = force * dx;
        *fy = force * dy;
        *fz = force * dz;
      } else {
        if (dx != 0)
          *fx = force / sqrt(1 + (dy * dy + dz * dz) / (dx * dx));
        else
          *fx = 0.0;
        if (dy != 0)
          *fy = force / sqrt(1 + (dx * dx + dz * dz) / (dy * dy));
        else
          *fy = 0.0;
        if (dz != 0)
          *fz = force / sqrt(1 + (dx * dx + dy * dy) / (dz * dz));
        else
          *fz = 0.0;
        //                if (dx < 0) *fx = -*fx;
        //                if (dy < 0) *fy = -*fy;
        //                if (dz < 0) *fz = -*fz;
      }
    }
    
    // ----------------------------------------------------------------------
    
    static inline double cutoff(double x){
      /*  cubic: return x * x * (2.0 * x - 3.0) + 1.0; */
      /*  quintic: return -6 * pow(x, 5) + 15 * pow(x, 4) - 10 * pow(x, 3) + 1; */
      /* Seventh order spline */
      //      return 20 * pow(x, 7) - 70 * pow(x, 6) + 84 * pow(x, 5) - 35 * pow(x, 4) + 1;
      return (((20 * x - 70) * x + 84) * x - 35) * x * x * x * x + 1;
    }
    
    // ----------------------------------------------------------------------
    
    static inline double dcutoff(double x){
      /*  cubic: return (6.0 * x * x - 6.0 * x); */
      /*  quintic: return -30 * pow(x, 4) + 60 * pow(x, 3) - 30 * pow(x, 2); */
    
      /* Seventh order spline */
      //      return 140 * pow(x, 6) - 420 * pow(x, 5) + 420 * pow(x, 4) - 140 * pow(x, 3);
      return (((140 * x - 420) * x + 420) * x - 140) * x * x * x;
    }
    
    }    // namespace LAMMPS_NS
    