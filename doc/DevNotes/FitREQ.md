
## General principle

We model hydrogen bonds by two methods
1. we decrease repulsive of non-covalent interatomic pairwise potential ( i.e. Pauli repulsion, modeled by $A/r_{ij}^{12}$ in Lennard-Jones or by $A \exp(-2 a r_{ij})$ in Morse and Buckingham.
   * This is done by multiplying amplitude of repulsion $A$ by some constant, $A_H = A(1 + H_{ij})$ where $H_{ij} = H_i H_j$
   * This correction is in effect only for atoms which have opposite Hbond correction charges  $H_i, H_j$, which means $H_{ij}<0$.
2. we add electron pairs as dummy atoms at certain distance (~0.5A) from host Hbond-acceptor (e.g. electronegative O,N atoms). We transfer some charge (0.1-0.6e) from host atom to electron pair. 
   * The electron pair is typically closer to Hbond-donor 
   * This is supposed to simulate the angular dependence of hydrogen bonds
   
## Only for selected atoms

* Fit REQ only for atoms which have non-zero Hb correction ( like `if(abs(Hbond)<1e-300) continue;` )
   * Or we may make a flag `bool* bHBs;   if(bHBs[ia]) continue`
   * Or we may make Hb atoms selection `int* iHbs;  for(int ia in iHbs){ ... };`
   * Or we may make new temp array which lists only hydrogen-bond corrected atoms
* For all of this we need do non-Hbonded reference
   * we compute reference Hbond energy first (using non-corrected LJQ or MorseQ) and store it in temporary array `atoms->userData->Emodel0`
   * to evaluate error we simply substract `dE = Eref - Emodel0 - Ecor` where Emodel0 is pre-calculated and Ecor is calculated just for selected Hbonded atoms.
 
## Electron pairs - short-range 
   * Electron pairs have certain charge on itself ( which is subtracted from host atom)
   * However they can have also some short range function (e.g. $\exp(-a r_{ij})$ ) which is attractive for all Hbond-donors (electron depleted Hydrogen atoms). 
      * we can check this by `Hb[ja]>0` 
