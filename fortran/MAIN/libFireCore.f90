
module firecore_options
    implicit none
    integer :: ioff_S=1, ioff_T=1, ioff_Vna=1, ioff_Vnl=1, ioff_Vxc=1, ioff_Vca=1, ioff_Vxc_ca=1, ioff_Ewald=1
    integer :: export_mode=0
end module firecore_options

! ==================================================
! ============ subroutine firecore_init
! ==================================================

subroutine firecore_setVerbosity( verbosity_, idebugWrite_ ) bind(c, name='firecore_setVerbosity')
    use iso_c_binding
    use options, only : verbosity, idebugWrite
    implicit none
    integer(c_int),intent(in), value :: verbosity_
    integer(c_int),intent(in), value :: idebugWrite_
    verbosity   = verbosity_
    idebugWrite = idebugWrite_
    write(*,*) "firecore_setVerbosity() ", verbosity, idebugWrite_
end subroutine

! see : https://stackoverflow.com/questions/29153501/when-i-pass-a-pointer-from-fortran-to-c-the-value-of-the-first-element-is-lost
subroutine firecore_getPointer_ratom( ratom_ ) bind(c, name='firecore_getPointer_ratom')
    use iso_c_binding
    use configuration
    implicit none
    type(C_PTR),      intent(out) :: ratom_
    ratom_ = c_loc( ratom ) 
end subroutine

subroutine firecore_getPointer_ftot( ftot_ ) bind(c, name='firecore_getPointer_ftot')
    use iso_c_binding
    use forces
    implicit none
    type(C_PTR),      intent(out) :: ftot_
    ftot_ = c_loc( ftot ) 
end subroutine 

subroutine firecore_getPointer_wfcoef( bbnkre_ ) bind(c, name='firecore_getPointer_wfcoef')
    use iso_c_binding
    use density
    implicit none
    type(C_PTR),      intent(out) :: bbnkre_
    bbnkre_ = c_loc( bbnkre ) 
end subroutine 

subroutine firecore_getPointer_charges( charges_ )  bind(c, name='firecore_getPointer_charges')
    use iso_c_binding
    use charges
    use options
    implicit none
    type(C_PTR),      intent(out) :: charges_
    if(verbosity.gt.0)write(*,*) "firecore_getPointer_charges() "
    if (iqout .eq. 2) then
        charges_ = c_loc( QMulliken_TOT )
    else
        charges_ = c_loc( QLowdin_TOT   )
    end if
end subroutine

subroutine firecore_set_lvs( lvs_ )  bind(c, name='firecore_set_lvs')
    use iso_c_binding
    use configuration
    use options
    implicit none
    real(c_double), dimension(3,3), intent(in) :: lvs_
    ! ========= body
    if(verbosity.gt.0)write(*,*) "firecore_set_lvs() "
    icluster = 0
    a1vec = lvs_(:,1)
    a2vec = lvs_(:,2)
    a3vec = lvs_(:,3)
end subroutine

subroutine firecore_initdir( )  bind(c, name='firecore_initdir' )
    use iso_c_binding
    use options
    use interactions
    use configuration
    use kpoints
    use loops
    use charges
    use integrals !, only : fdataLocation
    implicit none
    ! ========= body
    call initbasics()
    call readdata_mcweda ()
    call init_wfs(norbitals, nkpoints)
end subroutine

subroutine firecore_preinit( )  bind(c, name='firecore_preinit' )
    use iso_c_binding
    use options
    use configuration
    use kpoints
    use loops
    use charges
    use integrals !, only : fdataLocation
    use firecore_options
    implicit none
    ! ========= body
    !write(*,*) "DEBUG firecore_preinit() verbosity=", verbosity
    if(verbosity.gt.0)write(*,*) "firecore_preinit() "
    iparam_file = 0
    call initconstants ! (sigma, sigmaold, scf_achieved)
    scf_achieved = .true.
    call diagnostics (ioff2c, ioff3c, itestrange, testrange, ioff_S, ioff_T, ioff_Vna, ioff_Vnl, ioff_Vxc, ioff_Vca, ioff_Vxc_ca)    ! IF_DEF_DIAGNOSTICS
    !call readparam ()
    call set_default_params ()
    igrid = 1
    call checksum_options ()
    !write(*,*) "DEBUG firecore_preinit() END verbosity=", verbosity
    !write(*,*) "DEBUG firecore_preinit() icluster", icluster
end subroutine

subroutine firecore_init( natoms_, atomTypes, atomsPos ) bind(c, name='firecore_init')
    use iso_c_binding
    use options
    use loops
    use configuration
    use grid
    use interactions
    use integrals
    use density
    use kpoints
    use charges
    use FIRE, only: write_to_xyz
    implicit none
    ! ====== Parameters
    integer(c_int),                      intent(in), value  :: natoms_
    integer(c_int), dimension(natoms_),   intent(in)        :: atomTypes
    real(c_double), dimension(3,natoms_), intent(in)        :: atomsPos
    ! ====== global variables
    integer i, ispec, in1, iatom, numorbPP_max
    real distance
    real, dimension (3) :: vector
    logical zindata
    ! ====== Body
    !write(*,*) "DEBUG firecore_init() verbosity=", verbosity
    if(verbosity.gt.0)write(*,*) "firecore_init() "
    natoms = natoms_
    if(verbosity.gt.0) then
        write(*,*) "icluster", icluster
        write(*,*) "a1vec", a1vec
        write(*,*) "a3vec", a2vec
        write(*,*) "a2vec", a3vec
        write(*,*) "a2vec", a3vec
        write(*,*) "rm1,rm2,rm3,nrm ", rm1,rm2,rm3, nrm
        write(*,*) "em1,em2,em3,nem ", rm1,rm2,rm3, nem
    endif ! verbosity
 ! Allocate more arrays.
    allocate (degelec (natoms))
    allocate (iatyp (natoms))
    allocate (imass (natoms))
    allocate (ratom (3, natoms))
    allocate (nowMinusInitialPos (3, natoms))
    allocate (initialPosition (3, natoms))
    allocate (vatom (3, natoms))
    allocate (fragxyz(3, natoms))
    fragxyz(:,:) = 0
    allocate (symbol (natoms))
    allocate (xmass (natoms))
    allocate (ximage (3, natoms))
    ximage = 0.0d0
    iatyp(:) = atomTypes(:)
    ratom(:,:) = atomsPos(:,:)
    call readinfo ()   ! Read the info.dat file.  Allocate lsshPP, etc. inside!
    call cross (a2vec, a3vec, vector)
    Vouc=a1vec(1)*vector(1)+a1vec(2)*vector(2)+a1vec(3)*vector(3)
    call initboxes (1)
 ! Call make_munu. This routine determines all of the non-zero matrix elements for the two- and three-center matrix elements.  These non-zero matrix elements are determined based on selection rules.
    call make_munu   (nspecies)
    call make_munuPP (nspecies)
    call make_munuS  (nspecies)  
    call countOrbitals(numorb_max)                                                              ! IF_DEF_GRID_END


    call initcharges  (natoms, nspecies, itheory, ifixcharge, symbol) 
    call countIsorp()
    call countChargeDeglect(numorbPP_max)

    ! Initialise the grid variables.
    call initgrid ! (natoms, ia_g, ja_g, c2_g, rcut_g, ngrid)
    call get_info_orbital (natoms)
    if (itheory .eq. 2) call make_mu2shell (nspecies)
    call initamat(nspecies)
    allocate (xdot (0:5, 3, natoms)) ! Initialized below
 ! Allocate the stuff that depends on natoms, neigh_max, and numorb_max
    call allocate_neigh()
    call allocate_f() 
    call allocate_h() 
    call allocate_rho() 
    call readdata_mcweda ()
    call init_wfs(norbitals, nkpoints)

    !   usefull for projecting orbitals and density on grid
    call allocate_grid !(natoms, nspecies)
    call read_wf ()
    call read_vna ()

    call write_to_xyz( "#DEBUG libFireCore::firecore_init() ", 1 )
    return
end subroutine firecore_init

! ==================================================
! ============ subroutine firecore_SCF
! ==================================================

subroutine firecore_assembleH( iforce_, Kscf_, positions_ ) bind(c, name='firecore_assembleH')
    use iso_c_binding
    use options
    use configuration
    use debug
    use loops, only : Kscf
    implicit none
    ! ====== Parameters
    integer(c_int),                 intent(in),value :: iforce_
    integer(c_int),                 intent(in),value :: Kscf_
    real(c_double), dimension(3,natoms), intent(in) :: positions_
    ! ====== global variables
    integer ikpoint
    real, dimension (3) :: k_temp
    ! ====== Body
    if(verbosity.gt.0)write(*,*) "firecore_assembleH() "
    Kscf       = Kscf_
    iforce     = iforce_
    ratom(:,:) = positions_(:,:)
    call assemble_mcweda ()
    return
end subroutine

subroutine firecore_solveH( k_temp, ikpoint ) bind(c, name='firecore_solveH')
    use iso_c_binding
    use options
    use configuration
    use debug
    implicit none
    ! ====== Parameters
    integer(c_int),                 intent(in),value :: ikpoint
    real(c_double), dimension(3),   intent(in)       :: k_temp
    ! ====== global variables
    ! ====== Body
    if(verbosity.gt.0)write(*,*) "firecore_solveH() "
    call solveH ( ikpoint, k_temp )
    return
end subroutine

subroutine firecore_updateCharges( sigmatol_, sigma_ ) bind(c, name='firecore_updateCharges')
    use iso_c_binding
    use options
    use configuration
    use charges
    use loops, only: sigma,sigmatol
    use debug
    implicit none
    ! ====== Parameters
    real(c_double), intent(in) ,value   :: sigmatol_
    real(c_double), intent(out)         :: sigma_
    ! ====== global variables
    ! ====== Body
    if(verbosity.gt.0)write(*,*) "firecore_updateCharges() "
    sigmatol=sigmatol_
    call denmat ()
    sigma = sqrt(sum((Qin(:,:) - Qout(:,:))**2))
    call mixer ()
    sigma_=sigma
    return
end subroutine

subroutine firecore_SCF( nmax_scf, positions_, iforce_  )  bind(c, name='firecore_SCF')
    use iso_c_binding
    use options
    use loops
    use fire
    use configuration
    use energy
    use forces
    use interactions
    use integrals
    use density
    use kpoints
    use charges
    use debug
    implicit none
    ! ====== Parameters
    integer(c_int),                 intent(in),value :: nmax_scf
    integer(c_int),                 intent(in),value :: iforce_
    real(c_double), dimension(3,natoms), intent(in)  :: positions_
    ! ====== global variables
    integer ikpoint, i
    real, dimension (3) :: k_temp
    ! ====== Body
    if(verbosity.gt.0)write(*,*) "firecore_SCF() "
    ratom(:,:) = positions_(:,:)

    iforce  = iforce_
    ikpoint = 1
    scf_achieved = .false.
    max_scf_iterations = nmax_scf
    if(verbosity.gt.0)write(*,*) "!!!! SCF LOOP max_scf_iterations ", max_scf_iterations, scf_achieved
    do Kscf = 1, max_scf_iterations
        if(idebugWrite.gt.0)write(*,*) "! ======== Kscf ", Kscf
        call assemble_mcweda ()
        k_temp(:) = special_k(:,ikpoint)
        call solveH ( ikpoint, k_temp )
        call denmat ()
        sigma = sqrt(sum((Qin(:,:) - Qout(:,:))**2))
        if(verbosity.gt.0)write (*,*) "### SCF converged? ", scf_achieved, " Kscf ", Kscf, " |Qin-Oout| ",sigma," < tol ", sigmatol
        if( scf_achieved ) exit
        call mixer ()
    end do ! Kscf
    return
end subroutine firecore_SCF

subroutine firecore_evalForce( nmax_scf, positions_, forces_, energies, ixyzfile )  bind(c, name='firecore_evalForce')
    use iso_c_binding
    use options
    use loops
    use fire
    use configuration
    use energy
    use forces
    use interactions
    use integrals
    use density
    use kpoints
    use charges
    use debug
    implicit none
    ! ====== Parameters
    integer(c_int),                 intent(in),value :: nmax_scf
    real(c_double), dimension(3,natoms), intent(in)  :: positions_
    real(c_double), dimension(3,natoms), intent(out) :: forces_
    real(c_double), dimension(8),        intent(out) :: energies
    integer(c_int),                 intent(in),value :: ixyzfile
    ! ====== global variables
    integer ikpoint, i
    real, dimension (3) :: k_temp
    !real time_begin
    !real time_end
    ! ====== Body
    ! write(*,*) "DEBUG firecore_evalForce() verbosity=", verbosity
    if(verbosity.gt.0)write(*,*) "firecore_evalForce() "

    ratom(:,:) = positions_(:,:)
    iforce    = 1
    ftot(:,:) = 0
    ikpoint   = 1

    !do i=1,natoms  
    !    !write (*,*) "DEBUG firecore_evalForce() atom ", i, positions_(1,i), positions_(2,i), positions_(3,i)
    !    write (*,*) "DEBUG firecore_evalForce() atom ", i, iatyp(i), ratom(1,i), ratom(2,i), ratom(3,i)
    !end do

    scf_achieved = .false.
    max_scf_iterations = nmax_scf
    if(verbosity.gt.0)write(*,*) "!!!! SCF LOOP max_scf_iterations ", max_scf_iterations, scf_achieved
    do Kscf = 1, max_scf_iterations
        !if(idebugWrite.gt.0)write(*,*) "! ======== Kscf ", Kscf
        call assemble_mcweda ()
        !call debug_writeBlockedMat( "S_mat.log", s_mat )
        !call debug_writeBlockedMat( "H_mat.log", h_mat )
        k_temp(:) = special_k(:,ikpoint)
        call solveH ( ikpoint, k_temp )
        call denmat ()
        sigma = sqrt(sum((Qin(:,:) - Qout(:,:))**2))
        if(verbosity.gt.0)write (*,*) "### SCF converged? ", scf_achieved, " Kscf ", Kscf, " |Qin-Oout| ",sigma," < tol ", sigmatol
        if( scf_achieved ) exit
        !Qin(:,:) = Qin(:,:)*(1.0-bmix) + Qout(:,:)*bmix   ! linear mixer 
        call mixer ()
    end do ! Kscf
    !write(*,*) "DEBUG firecore_evalForce 1"
    call getenergy_mcweda ()
    !write(*,*) "DEBUG firecore_evalForce 2"
    call getforces_mcweda ()
    !write(*,*) "DEBUG firecore_evalForce 3"
    forces_(:,:) = ftot(:,:)
    !write(*,*) "DEBUG firecore_evalForce 4"
    energies(1) = etot    ! total energy
    energies(2) = ebs     ! band structure energy
    energies(3) = uiiuee     ! electrostatic
    energies(4) = etotxc_1c  ! exchange correlation on site
    energies(5) = uxcdcc     ! double counting correction ?
    energies(6) = atomic_energy
    energies(7) = efermi
    !write(*,*) "DEBUG firecore_evalForce 5"
    !call cpu_time (time_end)
    !write (*,*) ' FIREBALL RUNTIME : ',time_end-time_begin,'[sec]'
    !write(*,*) "DEBUG firecore_evalForce 6"
    if(ixyzfile .gt. 0) call write_to_xyz( "#firecore_evalForce() ", ixyzfile )
    !write(*,*) "DEBUG firecore_evalForce 7"
    if(verbosity.gt.0)write (*,*) '!!!! SCF LOOP DONE in ', Kscf, " iterations"
    return
end subroutine firecore_evalForce


subroutine firecore_relax( nstep, nmax_scf, positions_, forces_, fixPos, energies )  bind(c, name='firecore_relax')
    use iso_c_binding
    use options
    use loops
    use fire
    use configuration
    use energy
    use forces
    use interactions
    use integrals
    use density
    use kpoints
    use charges
    use grid
    use debug
    use timing
    implicit none
    ! ====== Parameters
    integer(c_int),                 intent(in),value :: nstep
    integer(c_int),                 intent(in),value :: nmax_scf
    real(c_double), dimension(3,natoms), intent(inout)  :: positions_
    real(c_double), dimension(3,natoms), intent(out) :: forces_
    integer(c_int), dimension(3,natoms), intent(out) :: fixPos
    real(c_double), dimension(8),        intent(out) :: energies
    !integer(c_int),                 intent(in),value :: ixyzfile
    ! ====== global variables
    integer i,j, in1
    integer ikpoint
    integer imu
    real, dimension (3) :: k_temp
    character (40) namewf
    ! ====== Body
    !write(*,*) "DEBUG firecore_relax() verbosity=", verbosity
    ikpoint = 1
    fragxyz(:,:) = fixPos(:,:)
    ratom(:,:)   = positions_(:,:) 
    forces_(:,:) = 0.0
    iparam_file = 1
    idebugWrite = 0
    iwrtxyz     = 1
    call init_FIRE( )
    call clean_ncall()
    do itime_step = 1,nstep
        iforce = 1
        scf_achieved = .false.
        do Kscf = 1, nmax_scf
            call assemble_mcweda()
            k_temp(:) = special_k(:,ikpoint)
            call solveH ( ikpoint, k_temp )
            call denmat ()
            sigma = sqrt(sum((Qin(:,:) - Qout(:,:))**2))
            if(verbosity.gt.1)write (*,*) "### SCF converged? ", scf_achieved, " Kscf ", Kscf, " |Qin-Oout| ",sigma," < tol ", sigmatol
            if( scf_achieved ) exit
            call mixer ()
        end do ! Kscf
        call getenergy_mcweda () 
        call getforces_mcweda ()
        call move_ions_FIRE (itime_step, iwrtxyz )   ! Move ions now
        call write_bas()
        !call write_to_xyz( "#firecore_evalForce() ", itime_step )
        if(verbosity.gt.0)write (*,'(A,i6,A,i6,A,f16.8,A,f16.8)') " ###### Time Step", itime_step,"(of ",nstepf,") |Fmax| ",  deltaFmax , " > force_tol " , force_tol
        if ( deltaFmax .lt. force_tol ) then
            if(verbosity.gt.0) write (*,*) "####### Relaxation DONE"
            exit
        endif
    end do ! istepf
    forces_   (:,:) = ftot(:,:)
    positions_(:,:) = ratom(:,:)
    energies(1) = etot    ! total energy
    energies(2) = ebs     ! band structure energy
    energies(3) = uiiuee     ! electrostatic
    energies(4) = etotxc_1c  ! exchange correlation on site
    energies(5) = uxcdcc     ! double counting correction ?
    energies(6) = atomic_energy
    energies(7) = efermi
end subroutine firecore_relax

subroutine firecore_getCharges( charges_ )  bind(c, name='firecore_getCharges')
    use iso_c_binding
    use configuration
    use charges
    use options
    real(c_double), dimension(natoms), intent(out) :: charges_
    if(verbosity.gt.0)write(*,*) "firecore_getCharges() "
    !write(*,*) "DEBUG firecore_getCharges START"
    if (iqout .eq. 2) then
        charges_(:) = QMulliken_TOT(:)
    else
        charges_(:) = QLowdin_TOT(:)
    end if
    !write(*,*) "DEBUG firecore_getCharges END"
end subroutine

subroutine firecore_get_Qin_shell( Qin_out ) bind(c, name='firecore_get_Qin_shell')
    use iso_c_binding
    use charges,       only: Qin
    use interactions,  only: nsh_max
    implicit none
    real(c_double), intent(out) :: Qin_out(*)
    integer iatom, issh, nsh, natoms
    nsh = size(Qin, 1)
    natoms = size(Qin, 2)
    do iatom = 1, natoms
        do issh = 1, nsh_max
            Qin_out((iatom-1)*nsh_max + issh) = 0.0d0
        end do
        do issh = 1, nsh
            Qin_out((iatom-1)*nsh_max + issh) = Qin(issh, iatom)
        end do
    end do
end subroutine

subroutine firecore_get_Qout_shell( Qout_out ) bind(c, name='firecore_get_Qout_shell')
    use iso_c_binding
    use charges,       only: Qout
    use interactions,  only: nsh_max
    implicit none
    real(c_double), intent(out) :: Qout_out(*)
    integer iatom, issh, nsh, natoms
    nsh = size(Qout, 1)
    natoms = size(Qout, 2)
    do iatom = 1, natoms
        do issh = 1, nsh_max
            Qout_out((iatom-1)*nsh_max + issh) = 0.0d0
        end do
        do issh = 1, nsh
            Qout_out((iatom-1)*nsh_max + issh) = Qout(issh, iatom)
        end do
    end do
end subroutine

! subroutine firecore_get_Qneutral_shell( Qneutral_out ) bind(c, name='firecore_get_Qneutral_shell')
!     use iso_c_binding
!     use charges,       only: Qneutral
!     implicit none
!     real(c_double), intent(out) :: Qneutral_out(*)
!     integer iatom, issh, nsh
!     nsh = size(Qneutral, 1)
!     do iatom = 1, size(Qneutral, 2)
!         Qneutral_out((iatom-1)*nsh + 1 : iatom*nsh) = 0.0d0
!         do issh = 1, nsh
!             Qneutral_out((iatom-1)*nsh + issh) = Qneutral(issh, iatom)
!         end do
!     end do
! end subroutine

subroutine firecore_get_Qneutral_shell( Qneutral_out ) bind(c, name='firecore_get_Qneutral_shell')
    use iso_c_binding
    use charges,       only: Qneutral, nzx
    use interactions,  only: nsh_max
    implicit none
    real(c_double), intent(out) :: Qneutral_out(*)
    integer ispec, issh, nsh, nspec, nspec_fdata
    
    ! Get actual dimensions from Fortran arrays
    nspec = 0
    nsh = 0
    if (allocated(Qneutral)) then
        nsh = size(Qneutral, 1)
        nspec = size(Qneutral, 2)
    end if
    if (allocated(nzx)) then
        nspec_fdata = size(nzx, 1)
    else
        nspec_fdata = nspec
    end if
    
    ! We don't know the exact size of Qneutral_out buffer passed from Python,
    ! but Python allocates it as dims.nsh_max * dims.nspecies_fdata.
    ! In OCL_Hamiltonian.py / FireCore.py, dims.nsh_max is usually 6 (matching interactions%nsh_max).
    ! dims.nspecies_fdata is usually size(nzx,1).
    
    ! Fill valid values from Qneutral. 
    ! We use nsh_max=6 as the stride to match Python's expectation of shape (6, nspecies_fdata)
    do ispec = 1, nspec_fdata
        do issh = 1, nsh_max
            Qneutral_out((ispec-1)*nsh_max + issh) = 0.0d0
        end do
    end do
    if (allocated(Qneutral)) then
        do ispec = 1, nspec
            do issh = 1, nsh
                Qneutral_out((ispec-1)*nsh_max + issh) = Qneutral(issh, ispec)
            end do
        end do
    end if
end subroutine

subroutine firecore_get_wfcoef( ikp, wfcoefs )  bind(c, name='firecore_get_wfcoef')
    use iso_c_binding
    use configuration
    use interactions
    use density
    use options
    integer(c_int)                                , intent(in),value :: ikp
    real(c_double), dimension(norbitals,norbitals), intent(out) :: wfcoefs
    if(verbosity.gt.0)write(*,*) "firecore_get_wfcoef() ", shape(bbnkre), ikp
    wfcoefs(:,:) = bbnkre(:,:,ikp)
end subroutine

subroutine firecore_set_wfcoef( iMO, ikp, wfcoefs )  bind(c, name='firecore_set_wfcoef')
    use iso_c_binding
    use configuration
    use interactions
    use density
    use options
    integer(c_int)                , intent(in),value :: ikp, iMO
    real(c_double), dimension(norbitals), intent(in) :: wfcoefs
    if(verbosity.gt.0)write(*,*) "firecore_set_wfcoef() ", shape(bbnkre), ikp, iMO
    !write (*,*) "shape(bbnkre) ", shape(bbnkre), ikp, iMO
    bbnkre(:,iMO,ikp) = wfcoefs(:)
end subroutine

subroutine firecore_setupGrid( Ecut_, ifixg0_, g0_, ngrid_, dCell  )  bind(c, name='firecore_setupGrid' )
    use iso_c_binding
    use grid
    use configuration
    use options
    implicit none
    ! ========= Parameters
    !integer(c_int),                 intent(in),value :: nmax_scf
    !real(c_double), dimension(3,natoms), intent(in) :: positions_
    !real(c_double), dimension(3,natoms), intent(out) :: forces_
    real   (c_double)               ,intent(in),value :: Ecut_
    integer(c_int)                  ,intent(in),value :: ifixg0_
    real   (c_double), dimension (3),intent(in) :: g0_
    integer(c_int),    dimension (3),intent  (inout) :: ngrid_
    real   (c_double), dimension (3,3),intent(inout) :: dCell
    !call readgrid !(iwrtewf)
    ! Namelist /mesh/ Ecut, iewform, npbands, pbands, ewfewin_max, ewfewin_min, ifixg0, g0
    ! ========= Variables
    logical bAutoGridSize
    ! ========= Body
    if(verbosity.gt.0)write(*,*) "firecore_setupGrid() "
    Ecut   = Ecut_
    ifixg0 = ifixg0_
    g0(:)  = g0_(:)
    call allocate_grid !(natoms, nspecies)
    call read_wf ()
    call read_vna ()
    bAutoGridSize = .False.
    if( (ngrid_(1).le.0) .or. (ngrid_(2).le.0) .or. (ngrid_(3).le.0) ) bAutoGridSize = .True.
    if( .not. bAutoGridSize ) then
        a1vec(:) = dCell(:,1)*ngrid_(1)
        a2vec(:) = dCell(:,2)*ngrid_(2)
        a3vec(:) = dCell(:,3)*ngrid_(3)
        ngrid(:)   = ngrid_(:)
    endif
    call initgrid_new( bAutoGridSize )
    ngrid_(1)=rm1
    ngrid_(2)=rm2
    ngrid_(3)=rm3
    dCell(:,:) = elvec(:,:)
    if(verbosity.gt.0)call write_grid_info()
end subroutine

subroutine firecore_getGridMO( iMO, ewfaux )  bind(c, name='firecore_getGridMO' )
    use iso_c_binding
    use grid
    use configuration
    !use options
    implicit none
    ! ========= Parameters
    integer(c_int), value :: iMO
    real(c_double), dimension (nrm), intent(out) :: ewfaux
    ! ========= Body
    !if(verbosity.gt.0)write(*,*) "firecore_getGridMO() ", iMO
    !allocate   ( ewfaux(0:nrm-1))
    !pewf => ewfaux
    write(*,*) "firecore_getGridMO ", iMO
    write(*,*) "DEBUG firecore_getGridMO 1 center_cell " 
    call center_cell ( .True. )
    write(*,*) "DEBUG firecore_getGridMO 2 project_orb iMO=", iMO 
    call project_orb( iMO, ewfaux )
end subroutine

subroutine firecore_getGridDens( ewfaux, f_den, f_den0 )  bind(c, name='firecore_getGridDens' )
    use iso_c_binding
    use grid
    use configuration
    implicit none
    ! ========= Parameters
    !integer(c_int), value :: imo0, imo1
    real(c_double), value :: f_den
    real(c_double), value :: f_den0
    real(c_double), dimension (nrm), intent(out) :: ewfaux
    ! ========= Body
    !if(verbosity.gt.0)write(*,*) "firecore_getGridDens() ", f_den, f_den0
    ewfaux = 0.0
    if( f_den*f_den > 1.e-16 ) then
        call project_dens( ewfaux, f_den )
        !call project_dens_new( ewfaux, f_den )
    end if 
    if( f_den0*f_den0 > 1.e-16 ) then
        call initdenmat0( )
        call project_dens0( ewfaux, f_den0 )
    end if 
end subroutine

subroutine firecore_orb2xsf( iMO )  bind(c, name='firecore_orb2xsf' )
    use iso_c_binding
    use options, only: icluster
    use grid
    implicit none
    ! ========= Parameters
    integer(c_int), value :: iMO
    ! ========= Variables
    real, dimension (:), allocatable :: ewfaux
    !real, target, dimension (:), allocatable :: ewfaux
    !real, dimension (:), pointer   :: pmat
    character(40)   :: namewf
    character(4)    :: name
    character (len=30) mssg
    integer i
    ! ========= Body
    !if(verbosity.gt.0)write(*,*) "firecore_orb2xsf() ", iMO
    allocate ( ewfaux(0:nrm-1))
    call project_orb(iMO,ewfaux)
    write (name,'(i4.4)') iMO
    namewf = 'bandplot_'//name//'.xsf'
    !pmat => ewfaux
    mssg = 'density_3D'
    !call writeout_xsf (namewf, mssg, pmat)
    call writeout_xsf (namewf, mssg, ewfaux)
    deallocate (ewfaux) 
end subroutine firecore_orb2xsf

subroutine firecore_dens2xsf( f_den0 )  bind(c, name='firecore_dens2xsf' )
    use iso_c_binding
    use options, only: icluster
    use grid
    implicit none
    ! ======= Local Parameters and Data Declaration
    real(c_double), value :: f_den0
    ! ========= Variables
    real, dimension (:), allocatable :: ewfaux
    character(40)   :: namewf
    character(4)    :: name
    character (len=30) mssg
    integer i
    ! ========= Body 
    !if(verbosity.gt.0)write(*,*) "firecore_dens2xsf() ", f_den0
    allocate ( ewfaux(0:nrm-1))
    ewfaux = 0.0
    call project_dens( ewfaux, 1.0 )
    if( f_den0*f_den0 > 1.e-16 ) then
        call initdenmat0( )
        call project_dens0( ewfaux, f_den0 )
    end if 
    namewf = 'FC_density.xsf'
    mssg = 'density_3D'
    call writeout_xsf (namewf, mssg, ewfaux)
    deallocate (ewfaux)
end subroutine firecore_dens2xsf

subroutine firecore_orb2points( iband,ikpoint, npoints, points, ewfaux )  bind(c, name='firecore_orb2points' )
    use iso_c_binding
    use options, only: icluster
    !use grid
    implicit none
    integer(c_int), value ::  iband, ikpoint, npoints
    real(c_double), dimension (3,npoints), intent (out) :: points
    real(c_double), dimension (  npoints), intent (out) :: ewfaux
    !if(verbosity.gt.0)write(*,*) "firecore_orb2points() "
    call project_orb_points( iband, ikpoint, npoints, points, ewfaux )
end subroutine firecore_orb2points

subroutine firecore_dens2points(npoints, points, f_den, f_den0, ewfaux_out) bind(c, name='firecore_dens2points')
    use iso_c_binding
    use configuration  ! For natoms, ratom, imass, num_orb, xl, neigh_*, etc.
    use density        ! For rho, rhoA
    use interactions   ! For Rc_max, nssh, lssh
    ! project_dens_points is in project_dens module
    !use project_dens, only: project_dens_points
    ! project_dens0_points and initdenmat0 are in project_dens0 module
    !use project_dens0, only: project_dens0_points, initdenmat0
    use options, only: verbosity
    implicit none

    integer(c_int), value, intent(in) :: npoints
    real(c_double), dimension(3, npoints), intent(in) :: points
    real(c_double), value, intent(in) :: f_den    ! Coefficient for rho_SCF
    real(c_double), value, intent(in) :: f_den0   ! Coefficient for rho_NA
    real(c_double), dimension(npoints), intent(out) :: ewfaux_out

    !if(verbosity.gt.0) 
    write(*,*) "firecore_dens2points() "
    if(verbosity.gt.0) write(*,*) "firecore_dens2points() npoints=", npoints, "f_den=", f_den, "f_den0=", f_den0

    ewfaux_out = 0.0_c_double

    if (npoints == 0) return

    if (abs(f_den) > 1.0e-16_c_double) then
        call project_dens_points(npoints, points, ewfaux_out, f_den)
    end if

    if (abs(f_den0) > 1.0e-16_c_double) then
        call initdenmat0() ! Ensure rhoA is initialized
        call project_dens0_points(npoints, points, ewfaux_out, f_den0)
    end if

end subroutine firecore_dens2points

subroutine firecore_getpsi( l, m, in1, issh, n, poss, ys )  bind(c, name='firecore_getpsi' )
    use iso_c_binding
    implicit none
    integer(c_int), value :: in1, issh, n, l, m
    !real(c_double), value :: x0, dx, theta, phi
    real(c_double), dimension (3,n), intent(out) :: poss
    real(c_double), dimension (n)  , intent(out) :: ys
    integer i
    real psi, dpsi !, x
    real, dimension (5)   :: Y      ! psi(x)
    real, dimension (3,5) :: dY     ! dpsi/dx
    !if(verbosity.gt.0)write(*,*) "firecore_getpsi() "
    do i=1,n
        call getpsi(in1,issh,NORM2(poss(:,i)),psi,dpsi)
        call getYlm(l,poss(:,i),Y,dY) 
        ys(i) = psi  *Y(m)
        !write (*,*) "", in1,issh,l,m, i,x,psi,Y(m)
    end do 
end subroutine

! ==================================================
! ============ Export H and S matrices
! ==================================================

subroutine firecore_get_HS_dims( &
    natoms_out, norbitals_out, nspecies_out, neigh_max_out, numorb_max_out, &
    nsh_max_out, ME2c_max_out, max_mu_dim1_out, max_mu_dim2_out, max_mu_dim3_out, mbeta_max_out, nspecies_fdata_out, nelec_out, neighPP_max_out &
  ) bind(c, name='firecore_get_HS_dims')
    use iso_c_binding
    use configuration, only: natoms, nspecies, xl ! nspecies here is for distinct species in current system
    use interactions,  only: norbitals, numorb_max, nsh_max, ME2c_max, mu
    use neighbor_map,  only: neigh_max, neighPP_max
    use charges,       only: nzx, ztot  ! nzx: species from Fdata; ztot: total electron count
    use options,       only: verbosity
    implicit none
    integer(c_int), intent(out) :: natoms_out, norbitals_out, nspecies_out, neigh_max_out, numorb_max_out
    integer(c_int), intent(out) :: nsh_max_out, ME2c_max_out
    integer(c_int), intent(out) :: max_mu_dim1_out, max_mu_dim2_out, max_mu_dim3_out, mbeta_max_out, nspecies_fdata_out, nelec_out, neighPP_max_out

    write (*,*) "firecore_get_HS_dims() verbosity=", verbosity

    natoms_out     = natoms
    ! Fortran side debug print for dimensions
! Fortran side debug print for dimensions
    if(verbosity.gt.1) then
        write(*,*) "firecore_get_HS_dims(): natoms     =", natoms
        write(*,*) "firecore_get_HS_dims(): norbitals  =", norbitals
        write(*,*) "firecore_get_HS_dims(): nspecies   =", nspecies
        write(*,*) "firecore_get_HS_dims(): neigh_max  =", neigh_max
        write(*,*) "firecore_get_HS_dims(): numorb_max =", numorb_max
        write(*,*) "firecore_get_HS_dims(): nsh_max    =", nsh_max
        if (allocated(mu)) then
            write(*,*) "firecore_get_HS_dims(): shape(mu)=", shape(mu)
        else
            if(verbosity.gt.0) write(*,*) "firecore_get_HS_dims(): xl not allocated yet!"
        end if
        if (allocated(nzx)) then
            write(*,*) "firecore_get_HS_dims(): shape(charges.nzx)=", shape(nzx)
        else
            if(verbosity.gt.0) write(*,*) "firecore_get_HS_dims(): charges.nzx not allocated yet!"
        end if
    end if
    norbitals_out  = norbitals
    nspecies_out   = nspecies ! Number of distinct species in the current calculation
    neigh_max_out  = neigh_max
    numorb_max_out = numorb_max
    nsh_max_out    = nsh_max
    ME2c_max_out   = ME2c_max ! Example, can add more ME_max if needed
    if (allocated(mu)) then
        max_mu_dim1_out = size(mu,1)
        max_mu_dim2_out = size(mu,2)
        max_mu_dim3_out = size(mu,3)
    else
        max_mu_dim1_out = 0 ! Default or error indicator
        max_mu_dim2_out = 0
        max_mu_dim3_out = 0
    end if
    if (allocated(xl)) then
        mbeta_max_out = size(xl,2) ! Get the actual second dimension size
    else
        mbeta_max_out = 0 ! Default or error indicator
    end if
    if (allocated(nzx)) then
        nspecies_fdata_out = size(nzx,1) ! Number of species defined in Fdata/info.dat
    else
        nspecies_fdata_out = 0
    end if

    neighPP_max_out = neighPP_max

    ! Export total electron count as integer approximation of ztot
    nelec_out = nint(ztot)
end subroutine firecore_get_HS_dims

subroutine firecore_get_HS_neighs( &
    num_orb_out, degelec_out, iatyp_out, &
    lssh_out, mu_out, nu_out, mvalue_out, nssh_out, nzx_out, &
    neighn_out, neigh_j_out, neigh_b_out, xl_out &
  ) bind(c, name='firecore_get_HS_neighs')
    use iso_c_binding
    use configuration, only: natoms, xl, nspecies
    use interactions,  only: num_orb, degelec, iatyp, lssh, mu, nu, mvalue, nssh, nsh_max
    use neighbor_map,  only: neighn, neigh_j, neigh_b, neigh_max
    use charges,       only: nzx
    use options, only : verbosity
    implicit none
    integer(c_int), dimension(nspecies), intent(out) :: num_orb_out
    integer(c_int), dimension(natoms),   intent(out) :: degelec_out
    integer(c_int), dimension(natoms),   intent(out) :: iatyp_out
    integer(c_int), dimension(nsh_max, nspecies), intent(out) :: lssh_out
    integer(c_int), dimension(size(mu,1), size(mu,2), size(mu,3)), intent(out) :: mu_out
    integer(c_int), dimension(size(nu,1), size(nu,2), size(nu,3)), intent(out) :: nu_out
    integer(c_int), dimension(size(mvalue,1), size(mvalue,2), size(mvalue,3)), intent(out) :: mvalue_out
    integer(c_int), dimension(nspecies), intent(out) :: nssh_out
    integer(c_int), dimension(size(nzx,1)), intent(out) :: nzx_out 
    integer(c_int), dimension(natoms),   intent(out) :: neighn_out
    integer(c_int), dimension(neigh_max, natoms), intent(out) :: neigh_j_out
    integer(c_int), dimension(neigh_max, natoms), intent(out) :: neigh_b_out
    real(c_double), dimension(3,size(xl,2)), intent(out) :: xl_out

    if(verbosity.gt.1) then
    !    write(*,"(A,I0,A,I0,A,I0,A,I0,A,I0,1X,I0,1X,I0,1X,I0)") "Fortran DIMS: h_mat_out(",numorb_max,",",numorb_max,",",neigh_max,",",natoms,") vs src h_mat ",shape(h_mat)
    !    write(*,"(A,I0,A,I0,A,I0,A,I0,A,I0,1X,I0,1X,I0,1X,I0)") "Fortran DIMS: s_mat_out(",numorb_max,",",numorb_max,",",neigh_max,",",natoms,") vs src s_mat ",shape(s_mat)
        write(*,"(A,I0,A,I0)")                                  "Fortran DIMS: num_orb_out(",nspecies,") vs src num_orb ",shape(num_orb) ! num_orb is 1D
        write(*,"(A,I0,A,I0)")                                  "Fortran DIMS: degelec_out(",natoms,") vs src degelec ",shape(degelec) ! degelec is 1D
        write(*,"(A,I0,A,I0)")                                  "Fortran DIMS: iatyp_out(",natoms,") vs src iatyp ",shape(iatyp)     ! iatyp is 1D
        write(*,"(A,I0,A,I0,A,I0,1X,I0)")                       "Fortran DIMS: lssh_out(",nsh_max,",",nspecies,") vs src lssh ",shape(lssh)
        write(*,"(A,I0,A,I0,A,I0,A,I0,1X,I0,1X,I0)")            "Fortran DIMS: mu_out(",size(mu,1),",",size(mu,2),",",size(mu,3),") vs src mu ",shape(mu)
        write(*,"(A,I0,A,I0,A,I0,A,I0,1X,I0,1X,I0)")            "Fortran DIMS: nu_out(",size(nu,1),",",size(nu,2),",",size(nu,3),") vs src nu ",shape(nu)
        write(*,"(A,I0,A,I0,A,I0,A,I0,1X,I0,1X,I0)")            "Fortran DIMS: mvalue_out(",size(mvalue,1),",",size(mvalue,2),",",size(mvalue,3),") vs src mvalue ",shape(mvalue)
        write(*,"(A,I0,A,I0)")                                  "Fortran DIMS: nssh_out(",nspecies,") vs src nssh ",shape(nssh)     ! nssh is 1D
        write(*,"(A,I0,A,I0)")                                  "Fortran DIMS: nzx_out(",size(nzx,1),") vs src charges.nzx ",shape(nzx) ! nzx is 1D
        write(*,"(A,I0,A,I0)")                                  "Fortran DIMS: neighn_out(",natoms,") vs src neighn ",shape(neighn)   ! neighn is 1D
        write(*,"(A,I0,A,I0,A,I0,1X,I0)")                       "Fortran DIMS: neigh_j_out(",neigh_max,",",natoms,") vs src neigh_j ",shape(neigh_j)
        write(*,"(A,I0,A,I0,A,I0,1X,I0)")                       "Fortran DIMS: neigh_b_out(",neigh_max,",",natoms,") vs src neigh_b ",shape(neigh_b)
        write(*,"(A,I0,A,I0,1X,I0)")                            "Fortran DIMS: xl_out(3,",size(xl,2),") vs src xl ",shape(xl)
    end if
    

    if(allocated(num_orb)) num_orb_out = num_orb
    if(allocated(degelec)) degelec_out = degelec
    if(allocated(iatyp))   iatyp_out   = iatyp
    if(allocated(lssh))    lssh_out    = lssh
    if(allocated(mu))      mu_out      = mu
    if(allocated(nu))      nu_out      = nu
    if(allocated(mvalue))  mvalue_out  = mvalue
    if(allocated(nssh))    nssh_out    = nssh
    if(allocated(nzx))     nzx_out     = nzx
    if(allocated(neighn))  neighn_out  = neighn
    if(allocated(neigh_j)) neigh_j_out = neigh_j
    if(allocated(neigh_b)) neigh_b_out = neigh_b
    if(allocated(xl))      xl_out      = xl
end subroutine firecore_get_HS_neighs


subroutine firecore_get_neigh_back( neigh_back_out ) bind(c, name='firecore_get_neigh_back')
    use iso_c_binding
    use configuration, only: natoms
    use neighbor_map,  only: neigh_max, neigh_back
    implicit none
    ! Export as (neigh_max,natoms) so Python C-order (natoms,neigh_max) can receive it without reordering.
    integer(c_int), dimension(neigh_max, natoms), intent(out) :: neigh_back_out
    neigh_back_out = 0
    if (allocated(neigh_back)) neigh_back_out = transpose(neigh_back)
end subroutine firecore_get_neigh_back

subroutine firecore_get_HS_neighsPP( neighPPn_out, neighPP_j_out, neighPP_b_out ) bind(c, name='firecore_get_HS_neighsPP')
    use iso_c_binding
    use configuration, only: natoms
    use neighbor_map,  only: neighPPn, neighPP_j, neighPP_b, neighPP_max
    implicit none
    integer(c_int), intent(out) :: neighPPn_out(*)   ! [natoms]
    integer(c_int), intent(out) :: neighPP_j_out(*)  ! [natoms*neighPP_max]
    integer(c_int), intent(out) :: neighPP_b_out(*)  ! [natoms*neighPP_max]
    integer :: iatom, ipp
    if(.not. allocated(neighPPn)) then
        write(*,*) 'Error: neighPPn not allocated in firecore_get_HS_neighsPP'
        stop
    end if
    if(.not. allocated(neighPP_j)) then
        write(*,*) 'Error: neighPP_j not allocated in firecore_get_HS_neighsPP'
        stop
    end if
    if(.not. allocated(neighPP_b)) then
        write(*,*) 'Error: neighPP_b not allocated in firecore_get_HS_neighsPP'
        stop
    end if
    do iatom = 1, natoms
        neighPPn_out(iatom) = neighPPn(iatom)
        do ipp = 1, neighPP_max
            neighPP_j_out((iatom-1)*neighPP_max + ipp) = neighPP_j(ipp, iatom)
            neighPP_b_out((iatom-1)*neighPP_max + ipp) = neighPP_b(ipp, iatom)
        end do
    end do
end subroutine firecore_get_HS_neighsPP

subroutine firecore_get_HS_sparse( h_mat_out, s_mat_out ) bind(c, name='firecore_get_HS_sparse')
    use iso_c_binding
    use configuration, only: natoms, ratom, xl
    use options,       only: itheory, itestrange, testrange, verbosity
    use firecore_options, only: ioff_S, ioff_T, ioff_Vna, ioff_Vnl, ioff_Vxc, ioff_Vca, ioff_Vxc_ca, ioff_Ewald, export_mode
    use interactions,  only: h_mat, s_mat, numorb_max, num_orb, imass, t_mat, vna, vnl, vxc, vxc_1c, vca, vxc_ca, ewaldlr, ewaldsr
    use neighbor_map,  only: neigh_max, neighPP_max, neighn, neigh_j, neigh_b, neighPPn, neighPP_j, neighPP_b
    implicit none
    real(c_double), dimension(numorb_max, numorb_max, neigh_max, natoms), intent(out) :: h_mat_out
    real(c_double), dimension(numorb_max, numorb_max, neigh_max, natoms), intent(out) :: s_mat_out
    integer :: iatom, ineigh, jatom, mbeta, in1, in2, nmu, nnu
    integer :: ineighPP, ineigh0
    real(c_double) :: dx, dy, dz, dist

    if (export_mode .eq. 0) then
        if(allocated(h_mat)) h_mat_out = h_mat
        if(allocated(s_mat)) s_mat_out = s_mat
        return
    end if

    h_mat_out = 0.0d0
    s_mat_out = 0.0d0
    if (ioff_S .ne. 0) then
        if(allocated(s_mat)) s_mat_out = s_mat
    end if

    if (.not. allocated(neighn)) then
        write(*,*) 'Error: neighn not allocated in firecore_get_HS_sparse export_mode=1'
        stop
    end if

    do iatom = 1, natoms
        in1 = imass(iatom)
        nmu = num_orb(in1)
        do ineigh = 1, neighn(iatom)
            jatom = neigh_j(ineigh, iatom)
            mbeta = neigh_b(ineigh, iatom)
            in2 = imass(jatom)
            nnu = num_orb(in2)

            if (ioff_T .ne. 0) then
                h_mat_out(:nmu,:nnu,ineigh,iatom) = h_mat_out(:nmu,:nnu,ineigh,iatom) + t_mat(:nmu,:nnu,ineigh,iatom)
            end if
            if (ioff_Vna .ne. 0) then
                h_mat_out(:nmu,:nnu,ineigh,iatom) = h_mat_out(:nmu,:nnu,ineigh,iatom) + vna(:nmu,:nnu,ineigh,iatom)
            end if

            if (itheory .eq. 3) then
                if (ioff_Vxc .ne. 0) then
                    h_mat_out(:nmu,:nnu,ineigh,iatom) = h_mat_out(:nmu,:nnu,ineigh,iatom) + vxc(:nmu,:nnu,ineigh,iatom)
                end if
                if (ioff_Vca .ne. 0) then
                    h_mat_out(:nmu,:nnu,ineigh,iatom) = h_mat_out(:nmu,:nnu,ineigh,iatom) + vca(:nmu,:nnu,ineigh,iatom)
                end if
            else
                if (ioff_Vxc .ne. 0) then
                    h_mat_out(:nmu,:nnu,ineigh,iatom) = h_mat_out(:nmu,:nnu,ineigh,iatom) + vxc(:nmu,:nnu,ineigh,iatom) + vxc_1c(:nmu,:nnu,ineigh,iatom)
                end if
                if (itheory .eq. 1 .or. itheory .eq. 2) then
                    if (ioff_Vca .ne. 0) then
                        h_mat_out(:nmu,:nnu,ineigh,iatom) = h_mat_out(:nmu,:nnu,ineigh,iatom) + vca(:nmu,:nnu,ineigh,iatom)
                        if (ioff_Ewald .ne. 0) then
                             h_mat_out(:nmu,:nnu,ineigh,iatom) = h_mat_out(:nmu,:nnu,ineigh,iatom) + ewaldlr(:nmu,:nnu,ineigh,iatom) - ewaldsr(:nmu,:nnu,ineigh,iatom)
                        end if
                    end if
                    if (ioff_Vxc_ca .ne. 0) then
                        h_mat_out(:nmu,:nnu,ineigh,iatom) = h_mat_out(:nmu,:nnu,ineigh,iatom) + vxc_ca(:nmu,:nnu,ineigh,iatom)
                    end if
                end if
            end if

            if (itestrange .eq. 0) then
                dx = ratom(1,iatom) - (xl(1,mbeta) + ratom(1,jatom))
                dy = ratom(2,iatom) - (xl(2,mbeta) + ratom(2,jatom))
                dz = ratom(3,iatom) - (xl(3,mbeta) + ratom(3,jatom))
                dist = sqrt(dx*dx + dy*dy + dz*dz)
                if (dist .gt. testrange) then
                    h_mat_out(:nmu,:nnu,ineigh,iatom) = 0.0d0
                    if (ioff_S .ne. 0) s_mat_out(:nmu,:nnu,ineigh,iatom) = 0.0d0
                end if
            end if
        end do
    end do

    ! Note: buildh.f90 does NOT add vnl into h_mat (it is kept separate and added later in k-space / solver).
    ! Therefore we only add vnl here when export_mode>=2, to provide an explicit "full-H" export.
    if (export_mode .ge. 2) then
        if (ioff_Vnl .ne. 0) then
            do iatom = 1, natoms
                in1 = imass(iatom)
                nmu = num_orb(in1)
                do ineighPP = 1, neighPPn(iatom)
                    jatom = neighPP_j(ineighPP, iatom)
                    mbeta = neighPP_b(ineighPP, iatom)
                    in2 = imass(jatom)
                    nnu = num_orb(in2)
                    ineigh0 = 0
                    do ineigh = 1, neighn(iatom)
                        if (neigh_j(ineigh,iatom) .eq. jatom .and. neigh_b(ineigh,iatom) .eq. mbeta) then
                            ineigh0 = ineigh
                            exit
                        end if
                    end do
                    if (ineigh0 .gt. 0) then
                        h_mat_out(:nmu,:nnu,ineigh0,iatom) = h_mat_out(:nmu,:nnu,ineigh0,iatom) + vnl(:nmu,:nnu,ineighPP,iatom)
                    else
                        if (verbosity .gt. 0) write(*,*) 'Warning: VNL neighbor not found in neigh list ', iatom, jatom, mbeta
                    end if
                end do
            end do
        end if
    end if
end subroutine firecore_get_HS_sparse

subroutine firecore_get_rho_sparse( rho_out ) bind(c, name='firecore_get_rho_sparse')
    use iso_c_binding
    use configuration, only: natoms
    use density,       only: rho
    use interactions,  only: numorb_max
    use neighbor_map,  only: neigh_max
    implicit none
    real(c_double), dimension(numorb_max, numorb_max, neigh_max, natoms), intent(out) :: rho_out
    call denmat ()
    if( .not. allocated(rho)) write(*,*) "Error: rho not allocated in firecore_get_rho_sparse" 
    rho_out = rho
end subroutine firecore_get_rho_sparse

subroutine firecore_get_rho_off_sparse( rho_off_out ) bind(c, name='firecore_get_rho_off_sparse')
    use iso_c_binding
    use configuration, only: natoms
    use density,       only: rho_off
    use interactions,  only: numorb_max
    use neighbor_map,  only: neigh_max
    implicit none
    real(c_double), dimension(numorb_max, numorb_max, neigh_max, natoms), intent(out) :: rho_off_out
    if( .not. allocated(rho_off)) write(*,*) "Error: rho_off not allocated in firecore_get_rho_off_sparse" 
    rho_off_out = rho_off
end subroutine firecore_get_rho_off_sparse


! ============================================================== 
! Universal export of 4D interaction arrays (sparse blocks) for debugging.
! DEBUG : TO EXPORT For checking /pyBall/FireballOCL/OCL_Hamiltonian.py
! code selects which array to export; out4d always has shape (numorb_max,numorb_max,neigh_max,natoms).
! codes:
!   1=vca   2=ewaldsr   3=ewaldlr   4=dip
!   5=vna   6=t_mat     7=vxc      8=vxc_ca
!   9=vnl (NOTE: uses neighPP_max internally; mapped into neigh_max when possible)
!  10=vca_ontopl   11=vca_ontopr   12=vca_atom
!  13=ewaldsr_2c_atom   14=ewaldsr_2c_ontop   15=ewaldsr_3c
subroutine firecore_export_interaction4D( code, out4d ) bind(c, name='firecore_export_interaction4D')
    use iso_c_binding
    use configuration, only: natoms
    use options,       only: verbosity
    use options,       only: idebugWrite
    use interactions,  only: numorb_max, vca, ewaldsr, ewaldlr, dip, vna, t_mat, vxc, vxc_ca, vnl, num_orb, imass
    ! --------------------------
    ! DEBUG : TO EXPORT For checking /pyBall/FireballOCL/OCL_Hamiltonian.py
    ! VCA component buffers live only in MODULES/debug.f90 and are allocated/filled
    ! only when idebugWrite>0 (assemble_ca_2c.f90 contains the gated hook).
    use debug,         only: dbg_vca_ontopl=>vca_ontopl, dbg_vca_ontopr=>vca_ontopr, dbg_vca_atom=>vca_atom, &
                             dbg_ewaldsr_2c_atom=>ewaldsr_2c_atom, dbg_ewaldsr_2c_ontop=>ewaldsr_2c_ontop, dbg_ewaldsr_3c=>ewaldsr_3c
    ! --------------------------
    use neighbor_map,  only: neigh_max, neighn, neigh_j, neigh_b, neighPPn, neighPP_j, neighPP_b
    implicit none
    integer(c_int), intent(in), value :: code
    real(c_double), dimension(numorb_max, numorb_max, neigh_max, natoms), intent(out) :: out4d
    integer :: iatom, ineigh, ineighPP, ineigh0, jatom, mbeta, in1, in2, nmu, nnu

    out4d = 0.0d0
    select case(code)
        case(1)
            if(allocated(vca)) out4d = vca
        case(2)
            if(allocated(ewaldsr)) out4d = ewaldsr
        case(3)
            if(allocated(ewaldlr)) out4d = ewaldlr
        case(4)
            if(allocated(dip)) out4d = dip
        case(5)
            if(allocated(vna)) out4d = vna
        case(6)
            if(allocated(t_mat)) out4d = t_mat
        case(7)
            if(allocated(vxc)) out4d = vxc
        case(8)
            if(allocated(vxc_ca)) out4d = vxc_ca
        case(10)
            if(idebugWrite .gt. 0) then
                if(allocated(dbg_vca_ontopl)) out4d = dbg_vca_ontopl
            end if
        case(11)
            if(idebugWrite .gt. 0) then
                if(allocated(dbg_vca_ontopr)) out4d = dbg_vca_ontopr
            end if
        case(12)
            if(idebugWrite .gt. 0) then
                if(allocated(dbg_vca_atom)) out4d = dbg_vca_atom
            end if
        case(13)
            if(idebugWrite .gt. 0) then
                if(allocated(dbg_ewaldsr_2c_atom)) out4d = dbg_ewaldsr_2c_atom
            end if
        case(14)
            if(idebugWrite .gt. 0) then
                if(allocated(dbg_ewaldsr_2c_ontop)) out4d = dbg_ewaldsr_2c_ontop
            end if
        case(15)
            if(idebugWrite .gt. 0) then
                if(allocated(dbg_ewaldsr_3c)) out4d = dbg_ewaldsr_3c
            end if
        case(9)
            ! vnl uses neighPP lists; map into neigh_max neighbor slots when possible
            if(.not. allocated(vnl)) return
            if(.not. allocated(neighn)) then
                if (verbosity .gt. 0) write(*,*) 'firecore_export_interaction4D(code=9): neighn not allocated'
                return
            end if
            do iatom = 1, natoms
                in1 = imass(iatom)
                nmu = num_orb(in1)
                do ineighPP = 1, neighPPn(iatom)
                    jatom = neighPP_j(ineighPP, iatom)
                    mbeta = neighPP_b(ineighPP, iatom)
                    in2 = imass(jatom)
                    nnu = num_orb(in2)
                    ineigh0 = 0
                    do ineigh = 1, neighn(iatom)
                        if (neigh_j(ineigh,iatom) .eq. jatom .and. neigh_b(ineigh,iatom) .eq. mbeta) then
                            ineigh0 = ineigh
                            exit
                        end if
                    end do
                    if (ineigh0 .gt. 0) then
                        out4d(:nmu,:nnu,ineigh0,iatom) = vnl(:nmu,:nnu,ineighPP,iatom)
                    else
                        if (verbosity .gt. 0) write(*,*) 'firecore_export_interaction4D(code=9): VNL neighbor not found ', iatom, jatom, mbeta
                    end if
                end do
            end do
        case default
            if (verbosity .gt. 0) write(*,*) 'firecore_export_interaction4D(): unknown code=', code
    end select
end subroutine firecore_export_interaction4D



subroutine firecore_set_options( ioff_S_, ioff_T_, ioff_Vna_, ioff_Vnl_, ioff_Vxc_, ioff_Vca_, ioff_Vxc_ca_, ioff_Ewald_ ) bind(c, name='firecore_set_options')
    use iso_c_binding
    use options
    use firecore_options
    implicit none
    integer(c_int), intent(in), value :: ioff_S_
    integer(c_int), intent(in), value :: ioff_T_
    integer(c_int), intent(in), value :: ioff_Vna_
    integer(c_int), intent(in), value :: ioff_Vnl_
    integer(c_int), intent(in), value :: ioff_Vxc_
    integer(c_int), intent(in), value :: ioff_Vca_
    integer(c_int), intent(in), value :: ioff_Vxc_ca_
    integer(c_int), intent(in), value :: ioff_Ewald_
    ioff_S    = ioff_S_
    ioff_T    = ioff_T_
    ioff_Vna  = ioff_Vna_
    ioff_Vnl  = ioff_Vnl_
    ioff_Vxc  = ioff_Vxc_
    ioff_Vca  = ioff_Vca_
    ioff_Vxc_ca = ioff_Vxc_ca_
    ioff_Ewald = ioff_Ewald_
    return
end subroutine

subroutine firecore_set_export_mode( export_mode_ ) bind(c, name='firecore_set_export_mode')
    use iso_c_binding
    use firecore_options
    implicit none
    integer(c_int), intent(in), value :: export_mode_
    export_mode = export_mode_
end subroutine

subroutine firecore_get_eigen( ikpoint, eigen_out ) bind(c, name='firecore_get_eigen')
    use iso_c_binding
    use interactions, only: norbitals
    use density,      only: eigen_k
    implicit none
    integer(c_int), value                :: ikpoint
    real(c_double), dimension(norbitals), intent(out) :: eigen_out

    eigen_out(:) = eigen_k(1:norbitals, ikpoint)
end subroutine firecore_get_eigen

subroutine firecore_get_HS_k(kpoint_vec, Hk_out, Sk_out) bind(c, name='firecore_get_HS_k')
    use iso_c_binding
    use interactions, only: norbitals
    ! Assuming ktransform is made available or its core logic is callable
    implicit none
    real(c_double),    dimension(3), intent(in) :: kpoint_vec
    complex(c_double), dimension(norbitals,norbitals), intent(out) :: Hk_out
    complex(c_double), dimension(norbitals,norbitals), intent(out) :: Sk_out
    ! Local copies for ktransform
    complex(c_double), dimension(norbitals,norbitals) :: local_Hk, local_Sk

    ! Call the k-transformation.
    ! ktransform itself uses h_mat, s_mat, vnl from interactions module
    ! and other arrays from configuration, neighbor_map etc.
    ! Ensure firecore_assembleH has been called.
    call ktransform(kpoint_vec, norbitals, local_Sk, local_Hk)

    Hk_out = local_Hk
    Sk_out = local_Sk

end subroutine firecore_get_HS_k

! Helper to get nspecies for dimensioning on C side (already have norbitals)
subroutine firecore_get_nspecies(nspecies_out) bind(c, name='firecore_get_nspecies')
    use iso_c_binding
    use configuration, only: nspecies
    implicit none
    integer(c_int), intent(out) :: nspecies_out
    nspecies_out = nspecies
end subroutine firecore_get_nspecies

subroutine firecore_get_nsh_max(nsh_max_out) bind(c, name='firecore_get_nsh_max')
    use iso_c_binding
    use interactions, only: nsh_max
    implicit none
    integer(c_int), intent(out) :: nsh_max_out
    nsh_max_out = nsh_max
end subroutine firecore_get_nsh_max

!=============================================================== 
! HORRIFIC BOILERPLATE RELATED TO DEBUG : TO EXPORT For checking /pyBall/FireballOCL/OCL_Hamiltonian.py
!===============================================================


! Add other simple dimension getters if needed, e.g., for ME2c_max
subroutine firecore_get_ME2c_max(ME2c_max_out) bind(c, name='firecore_get_ME2c_max')
    use iso_c_binding
    use interactions, only: ME2c_max
    implicit none
    integer(c_int), intent(out) :: ME2c_max_out
    ME2c_max_out = ME2c_max
end subroutine firecore_get_ME2c_max

subroutine firecore_scanHamPiece2c( interaction, isub, in1, in2, in3, dR, applyRotation, sx_out ) bind(c, name='firecore_scanHamPiece2c')
    use iso_c_binding
    use interactions, only : numorb_max
    implicit none
    integer(c_int), intent(in), value :: interaction, isub, in1, in2, in3
    real(c_double), dimension(3), intent(in) :: dR
    integer(c_int), intent(in), value :: applyRotation
    real(c_double), dimension(numorb_max, numorb_max), intent(out) :: sx_out
    
    real distance
    real(c_double), dimension(3,3) :: eps
    real(c_double), dimension(3,3,3) :: deps
    real(c_double), dimension(numorb_max, numorb_max) :: sx, spx
    integer iforce, i, j
    
    iforce = 0
    distance = sqrt(sum(dR**2))
    
    if( applyRotation .ne. 0 ) then
        call epsilon( (/0.0,0.0,0.0/), real(dR), eps )
    else
        eps = 0.0
        do i=1,3
            eps(i,i) = 1.0
        end do
    endif
    
    call doscentros (interaction, isub, iforce, in1, in2, in3, distance, eps, deps, sx, spx)
    sx_out(:,:) = sx(:,:)
end subroutine

subroutine firecore_scanHamPiece2c_batch( interaction, isub, in1, in2, in3, npoints, dRs, applyRotation, sx_out ) bind(c, name='firecore_scanHamPiece2c_batch')
    use iso_c_binding
    use interactions, only : numorb_max
    implicit none
    integer(c_int), intent(in), value :: interaction, isub, in1, in2, in3
    integer(c_int), intent(in), value :: npoints
    real(c_double), intent(in)  :: dRs(3, npoints)
    integer(c_int), intent(in), value :: applyRotation
    real(c_double), intent(out) :: sx_out(numorb_max, numorb_max, npoints)
    
    real(c_double) :: distance
    real(c_double), dimension(3,3) :: eps
    real(c_double), dimension(3,3,3) :: deps
    real(c_double), dimension(numorb_max, numorb_max) :: sx, spx
    real(c_double), dimension(3) :: dR
    integer :: iforce, i, ip
    
    do ip = 1, npoints
        iforce = 0
        dR(:) = dRs(:, ip)
        distance = sqrt(sum(dR**2))
        
        if( applyRotation .ne. 0 ) then
            call epsilon( (/0.0,0.0,0.0/), real(dR), eps )
        else
            eps = 0.0
            do i=1,3
                eps(i,i) = 1.0
            end do
        endif
        deps = 0.0
        call doscentros (interaction, isub, iforce, in1, in2, in3, distance, eps, deps, sx, spx)
        sx_out(:,:,ip) = sx(:,:)
    end do
end subroutine firecore_scanHamPiece2c_batch

subroutine firecore_scanHamPiece3c( interaction, isorp, in1, in2, indna, dRj, dRk, applyRotation, bcnax_out ) bind(c, name='firecore_scanHamPiece3c')
    use iso_c_binding
    use interactions, only : numorb_max
    use configuration, only : nspecies
    implicit none
    integer(c_int), intent(in), value :: interaction, isorp, in1, in2, indna
    real(c_double), dimension(3), intent(in) :: dRj, dRk
    integer(c_int), intent(in), value :: applyRotation
    real(c_double), dimension(numorb_max, numorb_max), intent(out) :: bcnax_out
    
    real(c_double), dimension(3) :: r1, r2, rna, r21, rnabc, rhat, sighat
    real x, y, cost
    real(c_double), dimension(3,3) :: eps
    real(c_double), dimension(numorb_max, numorb_max) :: bcnax
    integer i, j, maxtype
    
    ! Basis 1 at (0,0,0), Basis 2 at dRj, Potential at dRk
    rna = dRk
    r1  = (/0.0,0.0,0.0/)
    r2  = dRj
    
    r21 = r2 - r1
    y   = sqrt(sum(r21**2))
    rnabc = rna - (r1 + r2)*0.5
    x = sqrt(sum(rnabc**2))
    
    if (y .gt. 1.0d-05) then
        sighat = r21 / y
    else
        sighat = (/0.0, 0.0, 1.0/)
    endif
    
    if (x .gt. 1.0d-05) then
        rhat = rnabc / x
    else
        rhat = (/0.0, 0.0, 0.0/)
    endif
    
    cost = sum(sighat * rhat)
    
    if( applyRotation .ne. 0 ) then
        call epsilon ( real(rhat), real(sighat), eps)
    else
        eps = 0.0
        do i=1,3
            eps(i,i) = 1.0
        end do
    endif
    
    maxtype = nspecies
    call trescentros (interaction, isorp, maxtype, in1, in2, indna, x, y, cost, eps, bcnax, nspecies)
    bcnax_out(:,:) = bcnax(:,:)
end subroutine

! Raw 3c interpolation (no Legendre sum, no mvalue*sint, no recover_3c, no rotate)
subroutine firecore_scanHamPiece3c_raw( interaction, isorp, in1, in2, indna, dRj, dRk, bcnalist_out ) bind(c, name='firecore_scanHamPiece3c_raw')
    use iso_c_binding
    use interactions
    use integrals
    use timing, only : interaction_glob
    use options, only : verbosity
    implicit none
    integer(c_int), intent(in), value :: interaction, isorp, in1, in2, indna
    real(c_double), dimension(3), intent(in) :: dRj, dRk
    real(c_double), dimension(5, ME3c_max), intent(out) :: bcnalist_out
    
    real(c_double), dimension(3) :: r1, r2, rna, r21, rnabc
    real x, y
    integer iME, index, nx, ny
    real hx, hy, xxmax, yymax
    real Q_L, dQ_Ldx, dQ_Ldy
    external :: interpolate_2d
    integer, save :: debug_count = 0
    integer, parameter :: debug_limit = 64
    ! ensure interaction_glob is set for interpolate_2d counters
    interaction_glob = interaction
    
    bcnalist_out = 0.0d0
    
    ! geometry (same as firecore_scanHamPiece3c)
    rna = dRk
    r1  = (/0.0d0,0.0d0,0.0d0/)
    r2  = dRj
    r21 = r2 - r1
    y   = sqrt(sum(r21**2))
    rnabc = rna - (r1 + r2)*0.5d0
    x = sqrt(sum(rnabc**2))
    
    if ((interaction .ne. 1) .and. (interaction .ne. 3)) return
    
    index = icon3c(in1,in2,indna)
    if (interaction .eq. 1) then
        hx = hx_bcna(isorp,index)
        hy = hy_bcna(isorp,index)
        nx = numx3c_bcna(isorp,index)
        ny = numy3c_bcna(isorp,index)
        xxmax = x3cmax_bcna(isorp,index)
        yymax = y3cmax_bcna(isorp,index)
    else
        hx = hx_den3(isorp,index)
        hy = hy_den3(isorp,index)
        nx = numx3c_den3(isorp,index)
        ny = numy3c_den3(isorp,index)
        xxmax = x3cmax_den3(isorp,index)
        yymax = y3cmax_den3(isorp,index)
    end if
    if (verbosity >= 3 .and. debug_count < debug_limit) then
        debug_count = debug_count + 1
        write(*,'(A,I5,2X,2F12.6,2X,2F10.6,2X,2I5)') '[DBG-RAW] ip=', 1, x, y, hx, hy, nx, ny
    end if
    if (x .gt. xxmax .or. y .gt. yymax .or. x .lt. 0 .or. y .lt. 0) return
    
    do iME = 1, index_max3c(in1,in2)
        if (interaction .eq. 1) then
            call interpolate_2d (x, y, 0, nx, ny, hx, hy, bcna_01(1,1,iME,isorp,index), Q_L, dQ_Ldx, dQ_Ldy)
            bcnalist_out(1,iME) = Q_L
            call interpolate_2d (x, y, 0, nx, ny, hx, hy, bcna_02(1,1,iME,isorp,index), Q_L, dQ_Ldx, dQ_Ldy)
            bcnalist_out(2,iME) = Q_L
            call interpolate_2d (x, y, 0, nx, ny, hx, hy, bcna_03(1,1,iME,isorp,index), Q_L, dQ_Ldx, dQ_Ldy)
            bcnalist_out(3,iME) = Q_L
            call interpolate_2d (x, y, 0, nx, ny, hx, hy, bcna_04(1,1,iME,isorp,index), Q_L, dQ_Ldx, dQ_Ldy)
            bcnalist_out(4,iME) = Q_L
            call interpolate_2d (x, y, 0, nx, ny, hx, hy, bcna_05(1,1,iME,isorp,index), Q_L, dQ_Ldx, dQ_Ldy)
            bcnalist_out(5,iME) = Q_L
        else
            call interpolate_2d (x, y, 0, nx, ny, hx, hy, den3_01(1,1,iME,isorp,index), Q_L, dQ_Ldx, dQ_Ldy)
            bcnalist_out(1,iME) = Q_L
            call interpolate_2d (x, y, 0, nx, ny, hx, hy, den3_02(1,1,iME,isorp,index), Q_L, dQ_Ldx, dQ_Ldy)
            bcnalist_out(2,iME) = Q_L
            call interpolate_2d (x, y, 0, nx, ny, hx, hy, den3_03(1,1,iME,isorp,index), Q_L, dQ_Ldx, dQ_Ldy)
            bcnalist_out(3,iME) = Q_L
            call interpolate_2d (x, y, 0, nx, ny, hx, hy, den3_04(1,1,iME,isorp,index), Q_L, dQ_Ldx, dQ_Ldy)
            bcnalist_out(4,iME) = Q_L
            call interpolate_2d (x, y, 0, nx, ny, hx, hy, den3_05(1,1,iME,isorp,index), Q_L, dQ_Ldx, dQ_Ldy)
            bcnalist_out(5,iME) = Q_L
        end if
    end do
end subroutine firecore_scanHamPiece3c_raw

subroutine firecore_scanHamPiece3c_raw_batch( interaction, isorp, in1, in2, indna, npoints, dRjs, dRks, bcnalist_out ) bind(c, name='firecore_scanHamPiece3c_raw_batch')
    use iso_c_binding
    use interactions
    use integrals
    use timing, only : interaction_glob
    use options, only : verbosity
    implicit none
    integer(c_int), intent(in), value :: interaction, isorp, in1, in2, indna
    integer(c_int), intent(in), value :: npoints
    real(c_double), intent(in)  :: dRjs(3, npoints)
    real(c_double), intent(in)  :: dRks(3, npoints)
    real(c_double), intent(out) :: bcnalist_out(5, ME3c_max, npoints)
    
    real(c_double), dimension(3) :: r1, r2, rna, r21, rnabc
    real x, y
    integer iME, ip, index, nx, ny
    real hx, hy, xxmax, yymax
    real Q_L, dQ_Ldx, dQ_Ldy
    external :: interpolate_2d
    integer, save :: debug_count = 0
    integer, parameter :: debug_limit = 64
    
    bcnalist_out = 0.0d0
    
    if (interaction .ne. 1) return
    
    do ip = 1, npoints
        interaction_glob = interaction
        
        rna = dRks(:, ip)
        r1  = (/0.0d0,0.0d0,0.0d0/)
        r2  = dRjs(:, ip)
        r21 = r2 - r1
        y   = sqrt(sum(r21**2))
        rnabc = rna - (r1 + r2)*0.5d0
        x = sqrt(sum(rnabc**2))
        
        index = icon3c(in1,in2,indna)
        hx = hx_bcna(isorp,index)
        hy = hy_bcna(isorp,index)
        nx = numx3c_bcna(isorp,index)
        ny = numy3c_bcna(isorp,index)
        xxmax = x3cmax_bcna(isorp,index)
        yymax = y3cmax_bcna(isorp,index)
        if (verbosity >= 3 .and. debug_count < debug_limit) then
            debug_count = debug_count + 1
            write(*,'(A,I5,2X,2F12.6,2X,2F10.6,2X,2I5)') '[scanHamPiece3c_raw_batch] ip=', ip, x, y, hx, hy, nx, ny
        end if
        if (x .gt. xxmax .or. y .gt. yymax .or. x .lt. 0 .or. y .lt. 0) cycle
        
        do iME = 1, index_max3c(in1,in2)
            call interpolate_2d (x, y, 0, nx, ny, hx, hy, bcna_01(1,1,iME,isorp,index), Q_L, dQ_Ldx, dQ_Ldy)
            bcnalist_out(1,iME,ip) = Q_L
            call interpolate_2d (x, y, 0, nx, ny, hx, hy, bcna_02(1,1,iME,isorp,index), Q_L, dQ_Ldx, dQ_Ldy)
            bcnalist_out(2,iME,ip) = Q_L
            call interpolate_2d (x, y, 0, nx, ny, hx, hy, bcna_03(1,1,iME,isorp,index), Q_L, dQ_Ldx, dQ_Ldy)
            bcnalist_out(3,iME,ip) = Q_L
            call interpolate_2d (x, y, 0, nx, ny, hx, hy, bcna_04(1,1,iME,isorp,index), Q_L, dQ_Ldx, dQ_Ldy)
            bcnalist_out(4,iME,ip) = Q_L
            call interpolate_2d (x, y, 0, nx, ny, hx, hy, bcna_05(1,1,iME,isorp,index), Q_L, dQ_Ldx, dQ_Ldy)
            bcnalist_out(5,iME,ip) = Q_L
        end do
    end do
end subroutine firecore_scanHamPiece3c_raw_batch

! Export raw 3c table (bcna_0theta) for a given element pair / dn atom.
! Data are returned in Fortran (x-major, then y) order flattened into buffer.
subroutine firecore_export_bcna_table( interaction, isorp, in1, in2, indna, itheta, iME, maxsize, nx_out, ny_out, hx_out, hy_out, buffer, status ) bind(c, name='firecore_export_bcna_table')
    use iso_c_binding
    use interactions
    use integrals
    implicit none
    
    integer(c_int), intent(in), value :: interaction, isorp, in1, in2, indna, itheta, iME, maxsize
    integer(c_int), intent(out) :: nx_out, ny_out, status

    real(c_double), intent(out) :: hx_out, hy_out
    real(c_double), intent(out) :: buffer(maxsize)
    integer :: index, nx, ny, ix, iy, idx
    
    status = 0
    buffer = 0.0d0
    nx_out = 0; ny_out = 0
    hx_out = 0.0d0; hy_out = 0.0d0
    
    !if (interaction .ne. 1) then
    !    status = 1; return
    !end if
    
    index = icon3c(in1, in2, indna)
    nx = numx3c_bcna(isorp,index)
    ny = numy3c_bcna(isorp,index)
    hx_out = hx_bcna(isorp,index)
    hy_out = hy_bcna(isorp,index)
    
    if (iME < 1 .or. iME > index_max3c(in1,in2)) then
        status = 2; return
    end if
    if (itheta < 1 .or. itheta > 5) then
        status = 3; return
    end if
    if (maxsize < nx*ny) then
        status = 4; return
    end if
    
    nx_out = nx
    ny_out = ny
    

        do iy = 1, ny
            do ix = 1, nx
                idx = (iy-1)*nx + ix
                select case(itheta)
                case(1)
                    buffer(idx) = bcna_01(ix,iy,iME,isorp,index)
                case(2)
                    buffer(idx) = bcna_02(ix,iy,iME,isorp,index)
                case(3)
                    buffer(idx) = bcna_03(ix,iy,iME,isorp,index)
                case(4)
                    buffer(idx) = bcna_04(ix,iy,iME,isorp,index)
                case(5)
                    buffer(idx) = bcna_05(ix,iy,iME,isorp,index)
                end select
            end do
        end do

end subroutine firecore_export_bcna_table

subroutine firecore_scanHamPiece3c_batch( interaction, isorp, in1, in2, indna, npoints, dRjs, dRks, applyRotation, bcnax_out ) bind(c, name='firecore_scanHamPiece3c_batch')
    use iso_c_binding
    use interactions, only : numorb_max
    use configuration, only : nspecies
    implicit none
    integer(c_int), intent(in), value :: interaction, isorp, in1, in2, indna
    integer(c_int), intent(in), value :: npoints
    real(c_double), intent(in)  :: dRjs(3, npoints)
    real(c_double), intent(in)  :: dRks(3, npoints)
    integer(c_int), intent(in), value :: applyRotation
    real(c_double), intent(out) :: bcnax_out(numorb_max, numorb_max, npoints)
    
    real(c_double), dimension(3) :: r1, r2, rna, r21, rnabc, rhat, sighat
    real x, y, cost
    real(c_double), dimension(3,3) :: eps
    real(c_double), dimension(numorb_max, numorb_max) :: bcnax
    integer i, ip, maxtype
    
    maxtype = nspecies
    do ip = 1, npoints
        rna = dRks(:, ip)
        r1  = (/0.0,0.0,0.0/)
        r2  = dRjs(:, ip)
        
        r21 = r2 - r1
        y   = sqrt(sum(r21**2))
        rnabc = rna - (r1 + r2)*0.5
        x = sqrt(sum(rnabc**2))
        
        if (y .gt. 1.0d-05) then
            sighat = r21 / y
        else
            sighat = (/0.0, 0.0, 1.0/)
        endif
        
        if (x .gt. 1.0d-05) then
            rhat = rnabc / x
        else
            rhat = (/0.0, 0.0, 0.0/)
        endif
        
        cost = sum(sighat * rhat)
        
        if( applyRotation .ne. 0 ) then
            call epsilon ( real(rhat), real(sighat), eps)
        else
            eps = 0.0
            do i=1,3
                eps(i,i) = 1.0
            end do
        endif
        
        call trescentros (interaction, isorp, maxtype, in1, in2, indna, x, y, cost, eps, bcnax, nspecies)
        bcnax_out(:,:,ip) = bcnax(:,:)
    end do
end subroutine firecore_scanHamPiece3c_batch


! ============================================================== 
! EVEN MORE HORRIFIC BOILERPLATE RELATED TO DEBUG : TO EXPORT For checking /pyBall/FireballOCL/OCL_Hamiltonian.py
! ==============================================================

subroutine firecore_set_avg_rho_diag( enable, iatom, jatom, mbeta ) bind(c, name='firecore_set_avg_rho_diag')
    use iso_c_binding
    use debug, only: diag_avg_rho_enable, diag_avg_rho_iatom, diag_avg_rho_jatom_target, diag_avg_rho_mbeta_target, diag_avg_rho_ineigh
    implicit none
    integer(c_int), intent(in), value :: enable
    integer(c_int), intent(in), value :: iatom
    integer(c_int), intent(in), value :: jatom
    integer(c_int), intent(in), value :: mbeta
    diag_avg_rho_enable = (enable .ne. 0)
    diag_avg_rho_iatom = iatom
    diag_avg_rho_jatom_target = jatom
    diag_avg_rho_mbeta_target = mbeta
    diag_avg_rho_ineigh = -1
end subroutine firecore_set_avg_rho_diag

subroutine firecore_set_vca_diag( enable, iatom, jatom, mbeta, isorp ) bind(c, name='firecore_set_vca_diag')
    use iso_c_binding
    use debug, only: diag_vca_enable, diag_vca_iatom, diag_vca_jatom, diag_vca_mbeta, diag_vca_isorp
    implicit none
    integer(c_int), intent(in), value :: enable
    integer(c_int), intent(in), value :: iatom
    integer(c_int), intent(in), value :: jatom
    integer(c_int), intent(in), value :: mbeta
    integer(c_int), intent(in), value :: isorp
    diag_vca_enable = (enable .ne. 0)
    diag_vca_iatom  = iatom
    diag_vca_jatom  = jatom
    diag_vca_mbeta  = mbeta
    diag_vca_isorp  = isorp
end subroutine firecore_set_vca_diag

subroutine firecore_get_vca_diag_state( enable, iatom, jatom, mbeta, isorp ) bind(c, name='firecore_get_vca_diag_state')
    use iso_c_binding
    use debug, only: diag_vca_enable, diag_vca_iatom, diag_vca_jatom, diag_vca_mbeta, diag_vca_isorp
    implicit none
    integer(c_int), intent(out) :: enable
    integer(c_int), intent(out) :: iatom
    integer(c_int), intent(out) :: jatom
    integer(c_int), intent(out) :: mbeta
    integer(c_int), intent(out) :: isorp
    enable = 0
    if (diag_vca_enable) enable = 1
    iatom = diag_vca_iatom
    jatom = diag_vca_jatom
    mbeta = diag_vca_mbeta
    isorp = diag_vca_isorp
end subroutine firecore_get_vca_diag_state

subroutine firecore_get_avg_rho_diag_state( enable, iatom, jatom, mbeta ) bind(c, name='firecore_get_avg_rho_diag_state')
    use iso_c_binding
    use debug, only: diag_avg_rho_enable, diag_avg_rho_iatom, diag_avg_rho_jatom_target, diag_avg_rho_mbeta_target
    implicit none
    integer(c_int), intent(out) :: enable
    integer(c_int), intent(out) :: iatom
    integer(c_int), intent(out) :: jatom
    integer(c_int), intent(out) :: mbeta
    enable = 0
    if (diag_avg_rho_enable) enable = 1
    iatom = diag_avg_rho_iatom
    jatom = diag_avg_rho_jatom_target
    mbeta = diag_avg_rho_mbeta_target
end subroutine firecore_get_avg_rho_diag_state

subroutine firecore_get_avg_rho_diag_meta( iatom, ineigh, in1, in2, jatom, mbeta ) bind(c, name='firecore_get_avg_rho_diag_meta')
    use iso_c_binding
    use debug, only: diag_avg_rho_iatom, diag_avg_rho_ineigh, diag_avg_rho_in1, diag_avg_rho_in2, diag_avg_rho_jatom, diag_avg_rho_mbeta
    implicit none
    integer(c_int), intent(out) :: iatom, ineigh, in1, in2, jatom, mbeta
    iatom  = diag_avg_rho_iatom
    ineigh = diag_avg_rho_ineigh
    in1    = diag_avg_rho_in1
    in2    = diag_avg_rho_in2
    jatom  = diag_avg_rho_jatom
    mbeta  = diag_avg_rho_mbeta
end subroutine firecore_get_avg_rho_diag_meta

subroutine firecore_get_avg_rho_diag_eps2c( eps_out ) bind(c, name='firecore_get_avg_rho_diag_eps2c')
    use iso_c_binding
    use debug, only: diag_avg_rho_eps2c
    implicit none
    real(c_double), intent(out) :: eps_out(*)
    integer i,j
    do j=1,3
        do i=1,3
            eps_out((j-1)*3 + i) = diag_avg_rho_eps2c(i,j)
        end do
    end do
end subroutine firecore_get_avg_rho_diag_eps2c

subroutine firecore_get_avg_rho_diag_sm( sm_out ) bind(c, name='firecore_get_avg_rho_diag_sm')
    use iso_c_binding
    use debug, only: diag_avg_rho_sm
    implicit none
    real(c_double), intent(out) :: sm_out(*)
    integer i,j
    do j=1,6
        do i=1,6
            sm_out((j-1)*6 + i) = diag_avg_rho_sm(i,j)
        end do
    end do
end subroutine firecore_get_avg_rho_diag_sm

subroutine firecore_get_avg_rho_diag_rhom2c( m_out ) bind(c, name='firecore_get_avg_rho_diag_rhom2c')
    use iso_c_binding
    use debug, only: diag_avg_rho_rhom2c
    implicit none
    real(c_double), intent(out) :: m_out(*)
    integer i,j
    do j=1,6
        do i=1,6
            m_out((j-1)*6 + i) = diag_avg_rho_rhom2c(i,j)
        end do
    end do
end subroutine firecore_get_avg_rho_diag_rhom2c

subroutine firecore_get_avg_rho_diag_rhom3c( m_out ) bind(c, name='firecore_get_avg_rho_diag_rhom3c')
    use iso_c_binding
    use debug, only: diag_avg_rho_rhom3c
    implicit none
    real(c_double), intent(out) :: m_out(*)
    integer i,j
    do j=1,6
        do i=1,6
            m_out((j-1)*6 + i) = diag_avg_rho_rhom3c(i,j)
        end do
    end do
end subroutine firecore_get_avg_rho_diag_rhom3c

subroutine firecore_get_avg_rho_diag_rhooff_3c( m_out ) bind(c, name='firecore_get_avg_rho_diag_rhooff_3c')
    use iso_c_binding
    use debug, only: diag_avg_rho_rhooff_3c
    implicit none
    real(c_double), intent(out) :: m_out(*)
    integer i,j
    do j=1,8
        do i=1,8
            m_out((j-1)*8 + i) = diag_avg_rho_rhooff_3c(i,j)
        end do
    end do
end subroutine firecore_get_avg_rho_diag_rhooff_3c

subroutine firecore_get_avg_rho_diag_rhooff_final( m_out ) bind(c, name='firecore_get_avg_rho_diag_rhooff_final')
    use iso_c_binding
    use debug, only: diag_avg_rho_rhooff_final
    implicit none
    real(c_double), intent(out) :: m_out(*)
    integer i,j
    do j=1,8
        do i=1,8
            m_out((j-1)*8 + i) = diag_avg_rho_rhooff_final(i,j)
        end do
    end do
end subroutine firecore_get_avg_rho_diag_rhooff_final

! VXC debug export - bulk getter for all Vxc microtest data
subroutine firecore_get_vxc_diag_state( enable, iatom, ineigh ) bind(c, name='firecore_get_vxc_diag_state')
    use iso_c_binding
    use debug, only: diag_vxc_enable, diag_vxc_iatom, diag_vxc_ineigh
    implicit none
    integer(c_int), intent(out) :: enable
    integer(c_int), intent(out) :: iatom
    integer(c_int), intent(out) :: ineigh
    enable = diag_vxc_enable
    iatom = diag_vxc_iatom
    ineigh = diag_vxc_ineigh
end subroutine firecore_get_vxc_diag_state

subroutine firecore_get_vxc_diag_data( denmx_out, den1x_out, sx_out, dens_out, densij_out, &
                                       muxc_out, dmuxc_out, muxcij_out, dmuxcij_out, bcxcx_out, &
                                       arho_on_out, arhoi_on_out, rho_on_out, rhoi_on_out, &
                                       vxc_ca_out, vxc_1c_out ) bind(c, name='firecore_get_vxc_diag_data')
    use iso_c_binding
    use debug, only: dbg_vxc_denmx, dbg_vxc_den1x, dbg_vxc_sx, &
                      dbg_vxc_dens, dbg_vxc_densij, &
                      dbg_vxc_muxc, dbg_vxc_dmuxc, dbg_vxc_muxcij, dbg_vxc_dmuxcij, dbg_vxc_bcxcx, &
                      dbg_vxc_arho_on, dbg_vxc_arhoi_on, dbg_vxc_rho_on, dbg_vxc_rhoi_on, &
                      dbg_vxc_vxc_ca, dbg_vxc_vxc_1c
    implicit none
    real(c_double), dimension(4,4), intent(out) :: denmx_out
    real(c_double), dimension(4,4), intent(out) :: den1x_out
    real(c_double), dimension(4,4), intent(out) :: sx_out
    real(c_double), dimension(2,2), intent(out) :: dens_out
    real(c_double), dimension(2,2), intent(out) :: densij_out
    real(c_double), dimension(2,2), intent(out) :: muxc_out
    real(c_double), dimension(2,2), intent(out) :: dmuxc_out
    real(c_double), dimension(2,2), intent(out) :: muxcij_out
    real(c_double), dimension(2,2), intent(out) :: dmuxcij_out
    real(c_double), dimension(4,4), intent(out) :: bcxcx_out
    real(c_double), dimension(2,2), intent(out) :: arho_on_out
    real(c_double), dimension(2,2), intent(out) :: arhoi_on_out
    real(c_double), dimension(4,4), intent(out) :: rho_on_out
    real(c_double), dimension(4,4), intent(out) :: rhoi_on_out
    real(c_double), dimension(4,4), intent(out) :: vxc_ca_out
    real(c_double), dimension(4,4), intent(out) :: vxc_1c_out
    
    denmx_out = dbg_vxc_denmx
    den1x_out = dbg_vxc_den1x
    sx_out    = dbg_vxc_sx
    dens_out  = dbg_vxc_dens
    densij_out= dbg_vxc_densij
    muxc_out  = dbg_vxc_muxc
    dmuxc_out = dbg_vxc_dmuxc
    muxcij_out= dbg_vxc_muxcij
    dmuxcij_out= dbg_vxc_dmuxcij
    bcxcx_out = dbg_vxc_bcxcx
    arho_on_out = dbg_vxc_arho_on
    arhoi_on_out = dbg_vxc_arhoi_on
    rho_on_out = dbg_vxc_rho_on
    rhoi_on_out = dbg_vxc_rhoi_on
    vxc_ca_out = dbg_vxc_vxc_ca
    vxc_1c_out = dbg_vxc_vxc_1c
end subroutine firecore_get_vxc_diag_data