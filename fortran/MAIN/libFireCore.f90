
! ==================================================
! ============ subroutine firecore_init
! ==================================================

subroutine firecore_setVerbosity( verbosity_, idebugWrite_ ) bind(c, name='firecore_setVerbosity')
    use iso_c_binding
    use options
    implicit none
    integer(c_int),intent(in), value :: verbosity_
    integer(c_int),intent(in), value :: idebugWrite_
    verbosity   = verbosity_
    idebugWrite = idebugWrite_
    write(*,*) "DEBUG firecore_setVerbosity() ", verbosity
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

subroutine firecore_preinit( )  bind(c, name='firecore_preinit' )
    use iso_c_binding
    use options
    use configuration
    use kpoints
    use loops
    use charges
    use integrals !, only : fdataLocation
    implicit none
    ! ========= body
    write(*,*) "DEBUG firecore_preinit() verbosity=", verbosity
    if(verbosity.gt.0)write(*,*) "firecore_preinit() "
    iparam_file = 0
    call initconstants ! (sigma, sigmaold, scf_achieved)
    scf_achieved = .true.
    call diagnostics (ioff2c, ioff3c, itestrange, testrange)    ! IF_DEF_DIAGNOSTICS
    !call readparam ()
    call set_default_params ()
    igrid = 1
    call checksum_options ()
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
    write(*,*) "DEBUG firecore_init() verbosity=", verbosity
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
    !write(*,*) "DEBUG firecore_evalForce() verbosity=", verbosity
    if(verbosity.gt.0)write(*,*) "firecore_evalForce() "
    ratom(:,:) = positions_(:,:)
    iforce    = 1
    ftot(:,:) = 0
    ikpoint   = 1
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
    call getenergy_mcweda ()
    call getforces_mcweda ()
    forces_(:,:) = ftot(:,:)
    energies(1) = etot    ! total energy
    energies(2) = ebs     ! band structure energy
    energies(3) = uiiuee     ! electrostatic
    energies(4) = etotxc_1c  ! exchange correlation on site
    energies(5) = uxcdcc     ! double counting correction ?
    energies(6) = atomic_energy
    energies(7) = efermi
    !call cpu_time (time_end)
    !write (*,*) ' FIREBALL RUNTIME : ',time_end-time_begin,'[sec]'
    if(ixyzfile .gt. 0) call write_to_xyz( "#firecore_evalForce() ", ixyzfile )
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
    write(*,*) "DEBUG firecore_relax() verbosity=", verbosity
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
    !write(*,*) "firecore_getGridMO ", iMO
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
