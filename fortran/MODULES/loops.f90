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
integer  itime_step
real     force_tol
real     dt


contains
! ==================================
! ================= SUBROUTINES
! ==================================

subroutine init_wfs( norb, nkpt )
	use kpoints
	use density
    implicit none
    ! ==== parameters
    integer, intent(in) :: norb
    integer, intent(in) :: nkpt 
	!inquire(file="kpoints.kpts", exist=file_exists ) 
    ! === Body
    nkpoints = nkpt
    ! --- wavefunction coefs
    ! k_temp(:) = 0
    !if (iqout .ne. 2) then 
    !write(*,*) "init_wfs norb, nkpt ",  norb, nkpt, nkpoints
    allocate (blowre (norb, norb, nkpoints))
    allocate (blowim (norb, norb, nkpoints))
    blowre = 0.0d0
    blowim = 0.0d0
    !endif
    allocate (bbnkre (norb, norb, nkpoints))
    allocate (bbnkim (norb, norb, nkpoints))
    allocate (eigen_k (norb, nkpoints))
    !write(*,*) "shape(eigen_k)", shape(eigen_k) 
    eigen_k = 0.0d0
    bbnkim = 0.0d0
    bbnkre = 0.0d0

    ! --- Kpoints
    allocate (special_k(3,nkpoints))
    allocate (weight_k (nkpoints))
    special_k(:,:) = 0
    weight_k (:)   = 1

end subroutine init_wfs

end module loops
