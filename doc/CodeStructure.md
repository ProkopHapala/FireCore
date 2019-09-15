
## Electronic Problem Options

* Core ?
    * `itheory`
    * `itheory_xc`
* Optional
    * `iks`         Kohn-Sham grid
    * `iqmmm`       Quantum-Mechanics/Molecular-Mechanins
    * `iordern`     Linear-Scaling solver  (O(n))
    * `ibias`       Extarnal potential ramp (e.g. Bias Voltage)
* Possibly remove ?
    * `iforce`      Derivatives of Energy(Hamiltonian,Den-Mat,S-Mat) to obtain forces on atoms
* Not sure 
    * `idipole`    Dipole correction; ( seems to be connected with QMMM )
    * `igauss`     Calculate matrix elements (H-mat,S-mat) using Gaussian functions
    * `ivdw`       van-der-Waals correction
    * `iharmonic`  Harmonic correction to energy ? ( Zero energy vibrations ? ) 
* Minor / IO
    * `iwrtfpieces`  Write individual components of Hamiltonian (?) (For debugging?)

## NOTES:
 * `/ASSEMBLERS/getenergy.f90` contains nice overview of different electronics structure options (`itheory` etc.)


## COMMENTING MODULES

* `END_DEF_*`
* `IF_DEF_*`
* `IF_DEF_*_END`






