 module loops
! This module provide information about SCF cycle
! ===========================================================================

! SCF step
integer  max_scf_iterations
integer  Kscf
logical  scf_achieved

real tempfe   ! feri smearing

! Ionic  steps (MD )

integer  nstepi

 end module loops
