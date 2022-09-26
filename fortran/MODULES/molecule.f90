module molecule
    implicit none
    ! ========== Atoms

    ! TODO:
    ! * find where is allocated which array ?
    !   * we can use ALLOCATIONS to understand size and intent of each array    
    !       * use regular expression like  `allocate\s*\(\s*rho`    `allocate(.*)?rho`



    integer :: natoms
    real, dimension (:, :), allocatable :: ratom
    real, dimension (:, :), allocatable :: vatom             ! IF_DEF_move_END
    real, dimension (:, :), allocatable :: fatom             ! IF_DEF_move_END
    real, dimension (:), allocatable    :: smass             ! IF_DEF_move_END
    real, dimension (:), allocatable    :: xmass                 ! IF_DEF_move_END
    integer, dimension (:, :), allocatable :: fragxyz        ! IF_DEF_move_END
    real, dimension (3,3) :: avec

    ! ========== Electronic Structure
    integer norbitals       ! Total number of orbitals in the Hamiltonian and overlap matrices.

    ! ---- Charges
    real, dimension (:, :), allocatable :: Qs

    ! ---- Hamiltonian
    real, dimension (:, :, :, :), allocatable :: h_mat
    real, dimension (:, :, :, :), allocatable :: s_mat
    real, dimension (:, :, :, :), allocatable :: sm_mat

    ! ---- wavefunction 
    real, dimension (:, :),    allocatable, target :: eigen_k
    real, dimension (:, :, :), allocatable, target :: bbnkre
    real, dimension (:, :, :), allocatable, target :: bbnkim
    real, dimension (:, :, :), allocatable, target :: blowim
    real, dimension (:, :, :), allocatable, target :: blowre

    ! ---- density
    real, dimension (:, :, :, :), allocatable :: rho
    real, dimension (:, :, :, :), allocatable :: rhoPP

    ! ========== Quationable ?????????????
    real, dimension (:, :), allocatable    :: ximage                 ! IF_DEF_PBC_END
    real, dimension (:, :, :), allocatable :: xdot
    real, dimension (:,:), allocatable     :: xl




    ! ================== Neighbors

! Neighbor mapping and storage.
         integer, dimension (:, :), allocatable :: neigh_b
         integer, dimension (:, :), allocatable :: neigh_j
         integer, dimension (:),    allocatable :: neighn

! Neighbor (Pseudopotential) mapping and storage.
! nPP
         integer, dimension (:, :), allocatable :: nPP_b
         integer, dimension (:, :), allocatable :: nPP_j
         integer, dimension (:, :), allocatable :: nPP_map
         integer, dimension (:),    allocatable :: nPPn
         integer, dimension (:),    allocatable :: nPP_self
! nPPx
         integer, dimension (:, :), allocatable :: nPPx_b
         integer, dimension (:, :), allocatable :: nPPx_j 
         integer, dimension (:, :), allocatable :: nPPx_map
         integer, dimension (:, :), allocatable :: nPPx_point
         integer, dimension (:),    allocatable :: nPPxn
         integer, dimension (:),    allocatable :: nPPx_self
! neighPP
         integer, dimension (:, :), allocatable :: neighPP_b
         integer, dimension (:, :), allocatable :: neighPP_j
         integer, dimension (:),    allocatable :: neighPPn

! ! IF_DEF_vdW
! ! Neighbor mapping and storage - for van der Waals interactions.
!          integer, dimension (:, :), allocatable :: neigh_b_vdw
!          integer, dimension (:, :), allocatable :: neigh_j_vdw
!          integer, dimension (:),    allocatable :: neighn_vdw
! ! END_DEF_vdW

! Common neighbor mapping and storage.
         integer, dimension (:, :, :), allocatable :: neigh_comb 
         integer, dimension (:, :, :), allocatable :: neigh_comj
         integer, dimension (:, :, :), allocatable :: neigh_com_ng
         integer, dimension (:, :),    allocatable :: neigh_comm 
         integer, dimension (:),       allocatable :: neigh_comn 

! Common neighbor (Pseudopotential) mapping and storage.
! 3. party PP  common pairs
         integer, dimension (:, :, :), allocatable :: neighPP_comb 
         integer, dimension (:, :, :), allocatable :: neighPP_comj 
         integer, dimension (:, :),    allocatable :: neighPP_comm 
         integer, dimension (:),       allocatable :: neighPP_comn 

! Back neighbor mapping 
         integer, dimension (:,:), allocatable :: neigh_back

! neigh_self is the m value for the "self m"
         integer, dimension (:),  allocatable :: neigh_self
         integer, dimension (:),  allocatable :: neighPP_self
!         integer, dimension (:),  allocatable :: neigh_vdw_self   ! IF_DEF_vdW_END

! total neighbor mapping (mapping neigh and neighPP together)
         integer, dimension(:),   allocatable :: neighn_tot
         integer, dimension(:,:), allocatable :: neighj_tot
         integer, dimension(:,:), allocatable :: neighb_tot
         integer                              :: num_neig_maxtot



! ======================================================
! =========== Hamiltonian / interactions
! ======================================================

         integer, dimension (:, :),    allocatable :: index_max2c
         integer, dimension (:, :),    allocatable :: index_max3c
         integer, dimension (:, :),    allocatable :: index_max2cDipY
         integer, dimension (:, :),    allocatable :: index_max2cDipX
         integer, dimension (:, :),    allocatable :: lssh
         integer, dimension (:, :, :), allocatable :: mu
         integer, dimension (:, :, :), allocatable :: mvalue
         integer, dimension (:),       allocatable :: nssh
         integer                                      nssh_tot
         integer, dimension (:, :, :), allocatable :: nu
         integer, dimension (:),       allocatable :: num_orb
         integer, dimension (3)                    :: num_orb_sh  !only up to d orbitals
         integer, dimension (:, :, :), allocatable :: muDipY
         integer, dimension (:, :, :), allocatable :: nuDipY
         integer, dimension (:, :, :), allocatable :: muDipX
         integer, dimension (:, :, :), allocatable :: nuDipX

! To get information of the orbitals, re loaded in getinfo_orbital
         integer, dimension (:), allocatable :: getmssh
         integer, dimension (:), allocatable :: getlssh
         integer, dimension (:), allocatable :: getissh
         integer, dimension (:), allocatable :: getiatom
 
! Placement of matrix element in Hamiltonian and overlap matrices.  
! The placement is based on which atom number is accessed and the number
! of orbitals of the atom.
         integer, dimension (:), allocatable :: degelec

! Needed for extended hubbard interactions.
         integer, dimension (:, :), allocatable :: mu2shell
 
! These variables are specifically for the Kleinmann-Bylander pseudo-potentials
         integer, dimension (:, :),    allocatable :: index_maxPP
         integer, dimension (:, :),    allocatable :: lsshPP
         integer, dimension (:, :, :), allocatable :: muPP
         integer, dimension (:, :, :), allocatable :: nuPP
         integer, dimension (:),       allocatable :: nsshPP
         integer, dimension (:),       allocatable :: num_orbPP


! These variables are specifically for spherical density approximation 
! used in OLSXC method
         integer, dimension (:, :),    allocatable :: index_maxS
         integer, dimension (:, :, :), allocatable :: muS
         integer, dimension (:, :, :), allocatable :: nuS
         integer, dimension (:, :, :), allocatable :: mvalueS

! This array stores the species information for each atom.
         integer, dimension (:), allocatable :: iatyp
         integer, dimension (:), allocatable :: imass
 
         real, dimension (:, :), allocatable :: cl_PP

! These arrays store the interactions for the Hamiltonian matrix. 
! These will be dimensioned according to (mu, nu, natoms, neigh_max),
! where mu, nu are the orbitals for iatom and its neighbor.
         real, dimension (:, :, :, :), allocatable :: h_mat
         real, dimension (:, :, :, :), allocatable :: s_mat
         real, dimension (:, :, :, :), allocatable :: sm_mat

         real, dimension (:, :, :, :), allocatable :: sVNL
         real, dimension (:, :, :, :), allocatable :: t_mat
         real, dimension (:, :, :, :), allocatable :: vna
         real, dimension (:, :, :, :), allocatable :: vnl
         real, dimension (:, :, :, :), allocatable :: vnl2c
         real, dimension (:, :, :, :), allocatable :: vnl3c
         real, dimension (:, :, :, :), allocatable :: vxc
         real, dimension (:, :, :, :), allocatable :: vxc_1c

! These arrays store interactions which are needed for the DOGS 
! contributions to the Hamiltonian matrix. 
         real, dimension (:, :, :, :), allocatable :: dip
         real, dimension (:, :),       allocatable :: ewald
         real, dimension (:, :, :, :), allocatable :: ewaldlr
         real, dimension (:, :, :, :), allocatable :: ewaldsr
         real, dimension (:, :, :, :), allocatable :: vca
         real, dimension (:, :, :, :), allocatable :: vxc_ca
!         real, dimension (:, :, :, :), allocatable :: ewaldqmmm    ! IF_DEF_QMMM_END

! IF_DEF_dipole
! Dipole with XYZ components
        real, dimension (:, :, :),       allocatable :: dipcm
        real, dimension (:, :, :, :, :), allocatable :: dipc

! ======================================================
! =========== FORCES
! ======================================================

! Ewald and long-range forces
         real, dimension (:, :, :), allocatable :: dewald
         real, dimension (:, :),    allocatable :: fewald
         real, dimension (:, :),    allocatable :: flrew

!        real, dimension (:, :), allocatable :: flrew_qmmm     ! IF_DEF_QMMM_END
!        real, dimension (:, :), allocatable :: ftot_dftd3     ! IF_DEF_DFTD3_END
 
! Derivatives of interactions
         real, dimension (:, :, :, :, :),    allocatable :: dipp
         real, dimension (:, :, :, :),       allocatable :: dippcm
         real, dimension (:, :, :, :, :, :), allocatable :: dippc
         real, dimension (:, :, :, :, :),    allocatable :: sp_mat
         real, dimension (:, :, :, :, :),    allocatable :: spm_mat
         real, dimension (:, :, :, :, :),    allocatable :: spVNL
         real, dimension (:, :, :, :, :),    allocatable :: tp_mat

! Forces
         real, dimension (:, :), allocatable, target :: ftot
         real, dimension (:, :), allocatable :: ftotold
         real, dimension (:, :), allocatable :: ftotnew

         real, dimension (:, :), allocatable :: dusr
         real, dimension (:, :), allocatable :: dxcv
         real, dimension (:, :), allocatable :: fro
         real, dimension (:, :), allocatable :: ft

         real, dimension (:, :), allocatable :: deltaF

         real, dimension (:, :, :), allocatable :: fana
         real, dimension (:, :, :), allocatable :: fanl
         real, dimension (:, :, :), allocatable :: faxc
         real, dimension (:, :, :), allocatable :: fotna
         real, dimension (:, :, :), allocatable :: fotnl
         real, dimension (:, :, :), allocatable :: fotxc

         real, dimension (:, :), allocatable :: f3naa, f3nab, f3nac
         real, dimension (:, :), allocatable :: f3nla, f3nlb, f3nlc
         real, dimension (:, :), allocatable :: f3xca, f3xcb, f3xcc

! Forces - OSLXC
         real, dimension (:, :, :), allocatable :: dxcdcc
         real, dimension (:, :, :, :, :), allocatable :: arhop_off
         real, dimension (:, :, :, :, :), allocatable :: arhopij_off 
         real, dimension (:, :, :, :, :), allocatable :: rhop_off
         real, dimension (:, :, :, :, :), allocatable :: rhopij_off
         real, dimension (:, :, :, :, :), allocatable :: rhop_on
         real, dimension (:, :, :, :, :), allocatable :: arhop_on

! Forces - DOGS
         real, dimension (:, :, :), allocatable :: faca
         real, dimension (:, :, :), allocatable :: faxc_ca
         real, dimension (:, :, :), allocatable :: fotca
         real, dimension (:, :, :), allocatable :: fotxc_ca

         real, dimension (:, :), allocatable :: f3caa, f3cab, f3cac
         real, dimension (:, :), allocatable :: f3xca_ca, f3xcb_ca, f3xcc_ca








end module
! ! END_DEF_QMMM
