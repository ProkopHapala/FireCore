module conf
  real(8), parameter :: tol = 0.05d0
  integer, save :: natoms, nbonds, nangles, ntorsions, ninversions
  integer, save :: natomtypes, nbondtypes, nangletypes, ntorsiontypes, ninversiontypes
  real(8), save :: side(3,3), side_lammps(3,3)
  integer, allocatable, save :: neighbors(:), connectivity(:,:), bonds(:,:), angles(:,:), torsions(:,:), inversions(:,:)
  integer, allocatable, save :: atomtypesmap(:), bondtypesmap(:), angletypesmap(:), torsiontypesmap(:), inversiontypesmap(:)
  integer, allocatable, save :: valence_target(:), bond_orders_int(:)
  real(8), allocatable, save :: pos(:,:), pos_lammps(:,:), q(:), bond_orders(:)
  logical, allocatable, save :: set_atom(:), set_bond(:)!, adhoc_atom(:), adhoc_bond(:)
  character(2), allocatable, save :: tipo(:)
  character(5), allocatable, save :: ufftypes(:)
end module conf

module constants
  real(8), parameter :: radius_H = 1.487d0 ! Ang
  real(8), parameter :: radius_C = 1.908d0 ! Ang
  real(8), parameter :: radius_N = 1.780d0 ! Ang
  real(8), parameter :: radius_O = 1.661d0 ! Ang
  real(8), parameter :: mass_H = 1.00794d0 ! amu
  real(8), parameter :: mass_C = 12.0107d0 ! amu
  real(8), parameter :: mass_N = 14.0067d0 ! amu
  real(8), parameter :: mass_O = 15.9994d0 ! amu
end module constants

module uff
  integer, save :: nufftypes
  real(8), allocatable, save :: r_I(:), theta_0(:), x_I(:), D_I(:), zeta(:), Zstar_I(:), chi(:), J_I(:), R(:), eta(:)
  real(8), allocatable, save :: V_I(:), U_I(:), chi_reax(:), eta_reax(:), gamma_reax(:)
  character(5), allocatable, save :: sym(:)
end module uff

module pars
  real(8), allocatable, save :: vdw_pars(:,:), bond_pars(:,:), angle_pars(:,:), torsion_pars(:,:), inversion_pars(:,:)
  character(15), allocatable, save :: angle_tags(:)
end module pars
