 module loops
! This module provide information about SCF cycle
! ===========================================================================

! SCF step
integer  max_scf_iterations
integer  Kscf
logical  scf_achieved

real tempfe   ! feri smearing

real bmix     ! charge mixing coeficient
real  sigmatol
real  sigmaold
real  sigma

!integer ialgmix

! idmix: Defines the mixing procedure: idmix=1 means simple mixing
! For larger idmix values, the choice of bmix becomes less important
 integer, parameter ::  idmix = 6

! Ionic  steps (MD )

integer  nstepf

 end module loops
