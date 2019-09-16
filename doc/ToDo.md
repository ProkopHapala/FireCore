
# To Do first

 * Define one option of `iTheory`, `idogs` etc., which functions in `ASSEMBLE` a `INTERACTIONS` are always essential to asseble hamiltonian (and overlap matrix)?
 * Try build something which assemble some minimalistic Hamiltonian

## Minor

 * Automatic script to reformat line breaks like this:
  ```
   rhoij_off(imu,inu,ineigh,iatom) = rhoij_off(imu,inu,ineigh,iatom) &
     &         + rhomx(imu,inu)*Qneutral(isorp,in2)
  ```
 * rename spanish names of files and subroutines `unocentros`, `doscentros`, `trescentros`

## Contemplation

 * Possibly use [**fypp**](https://github.com/aradi/fypp)