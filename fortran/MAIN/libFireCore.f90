
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
    implicit none
    ! ========= body
    !write(*,*) "DEBUG firecore_preinit() verbosity=", verbosity
    if(verbosity.gt.0)write(*,*) "firecore_preinit() "
    iparam_file = 0
    call initconstants ! (sigma, sigmaold, scf_achieved)
    scf_achieved = .true.
    call diagnostics (ioff2c, ioff3c, itestrange, testrange)    ! IF_DEF_DIAGNOSTICS
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
    nsh_max_out, ME2c_max_out, max_mu_dim1_out, max_mu_dim2_out, max_mu_dim3_out, mbeta_max_out, nspecies_fdata_out &
  ) bind(c, name='firecore_get_HS_dims')
    use iso_c_binding
    use configuration, only: natoms, nspecies, xl ! nspecies here is for distinct species in current system
    use interactions,  only: norbitals, numorb_max, nsh_max, ME2c_max, mu
    use neighbor_map,  only: neigh_max
    use charges, only: nzx  ! This nzx is dimensioned by total species in info.dat
    use options, only: verbosity
    implicit none
    integer(c_int), intent(out) :: natoms_out, norbitals_out, nspecies_out, neigh_max_out, numorb_max_out
    integer(c_int), intent(out) :: nsh_max_out, ME2c_max_out
    integer(c_int), intent(out) :: max_mu_dim1_out, max_mu_dim2_out, max_mu_dim3_out, mbeta_max_out, nspecies_fdata_out

    write (*,*) "firecore_get_HS_dims() verbosity=", verbosity

    natoms_out     = natoms
    ! Fortran side debug print for dimensions
! Fortran side debug print for dimensions
    if(verbosity.gt.1) then
        write(*,*) "Fortran firecore_get_HS_dims: natoms (current system)=", natoms
        write(*,*) "Fortran firecore_get_HS_dims: norbitals=", norbitals
        write(*,*) "Fortran firecore_get_HS_dims: nspecies (distinct in current system)=", nspecies
        write(*,*) "Fortran firecore_get_HS_dims: neigh_max=", neigh_max
        write(*,*) "Fortran firecore_get_HS_dims: numorb_max=", numorb_max
        write(*,*) "Fortran firecore_get_HS_dims: nsh_max=", nsh_max
        if (allocated(mu)) then
            write(*,*) "Fortran firecore_get_HS_dims: shape(mu)=", shape(mu)
        else
              if(verbosity.gt.0) write(*,*) "Fortran firecore_get_HS_dims: xl not allocated yet!"
        end if
        if (allocated(nzx)) then
            write(*,*) "Fortran firecore_get_HS_dims: shape(charges.nzx)=", shape(nzx)
        else
            if(verbosity.gt.0) write(*,*) "Fortran firecore_get_HS_dims: charges.nzx not allocated yet!"
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
end subroutine firecore_get_HS_dims

subroutine firecore_get_HS_sparse( &
    h_mat_out, s_mat_out, &
    num_orb_out, degelec_out, iatyp_out, &
    lssh_out, mu_out, nu_out, mvalue_out, nssh_out, nzx_out, &
    neighn_out, neigh_j_out, neigh_b_out, xl_out &
  ) bind(c, name='firecore_get_HS_sparse')
    use iso_c_binding
    use configuration, only: natoms, xl, nspecies ! nspecies is from configuration
    use interactions,  only: h_mat, s_mat, num_orb, degelec, iatyp, lssh, mu, nu, mvalue, nssh, numorb_max, nsh_max
    use neighbor_map,  only: neighn, neigh_j, neigh_b, neigh_max
    use charges, only: nzx
    use options, only: verbosity
    implicit none
    ! Output arrays (pointers to pre-allocated memory from C/Python)
    real(c_double), dimension(numorb_max, numorb_max, neigh_max, natoms), intent(out) :: h_mat_out
    real(c_double), dimension(numorb_max, numorb_max, neigh_max, natoms), intent(out) :: s_mat_out
    integer(c_int), dimension(nspecies), intent(out) :: num_orb_out
    integer(c_int), dimension(natoms),   intent(out) :: degelec_out
    integer(c_int), dimension(natoms),   intent(out) :: iatyp_out
    integer(c_int), dimension(nsh_max, nspecies), intent(out) :: lssh_out
    integer(c_int), dimension(size(mu,1), size(mu,2), size(mu,3)), intent(out) :: mu_out
    integer(c_int), dimension(size(nu,1), size(nu,2), size(nu,3)), intent(out) :: nu_out
    integer(c_int), dimension(size(mvalue,1), size(mvalue,2), size(mvalue,3)), intent(out) :: mvalue_out
    integer(c_int), dimension(nspecies), intent(out) :: nssh_out
    integer(c_int), dimension(size(nzx,1)), intent(out) :: nzx_out ! Use actual size of charges.nzx
    integer(c_int), dimension(natoms),   intent(out) :: neighn_out
    integer(c_int), dimension(neigh_max, natoms), intent(out) :: neigh_j_out
    integer(c_int), dimension(neigh_max, natoms), intent(out) :: neigh_b_out
    real(c_double), dimension(3,size(xl,2)), intent(out) :: xl_out

    write (*,*) "firecore_get_HS_sparse() verbosity=", verbosity," h_mat, s_mat allocated? ", allocated(h_mat), allocated(s_mat)

    ! Check if source arrays are allocated
    if (.not. allocated(h_mat)) then
        if(verbosity.gt.0) write(*,*) "firecore_get_HS_sparse: h_mat not allocated!"
        return
    end if
    if (.not. allocated(s_mat)) then
        if(verbosity.gt.0) write(*,*) "firecore_get_HS_sparse: s_mat not allocated!"
        return
    end if
    ! Add more checks for other arrays as needed

    if(verbosity.gt.1) then
        write(*,"(A,I4,A,I4,A,I4,A,I4,A,4I4)") "Fortran DIMS: h_mat_out(",numorb_max,",",numorb_max,",",neigh_max,",",natoms,") vs src h_mat",shape(h_mat)
        write(*,"(A,I4,A,I4,A,I4,A,I4,A,4I4)") "Fortran DIMS: s_mat_out(",numorb_max,",",numorb_max,",",neigh_max,",",natoms,") vs src s_mat",shape(s_mat)
        write(*,"(A,I2,A,I2)")       "Fortran DIMS: num_orb_out(",nspecies,") vs src num_orb",shape(num_orb)
        write(*,"(A,I2,A,I2)")       "Fortran DIMS: degelec_out(",natoms,") vs src degelec",shape(degelec)
        write(*,"(A,I2,A,I2)")       "Fortran DIMS: iatyp_out(",natoms,") vs src iatyp",shape(iatyp)
        write(*,"(A,I2,A,I2,A,I2)")  "Fortran DIMS: lssh_out(",nsh_max,",",nspecies,") vs src lssh",shape(lssh)
        write(*,"(A,I2,A,I2,A,I2,A,I2,A,I2,A,I2)") "Fortran DIMS: mu_out(",size(mu,1),",",size(mu,2),",",size(mu,3),") vs src mu",shape(mu)
        write(*,"(A,I2,A,I2,A,I2,A,I2,A,I2,A,I2)") "Fortran DIMS: nu_out(",size(nu,1),",",size(nu,2),",",size(nu,3),") vs src nu",shape(nu)
        write(*,"(A,I2,A,I2,A,I2,A,I2,A,I2,A,I2)") "Fortran DIMS: mvalue_out(",size(mvalue,1),",",size(mvalue,2),",",size(mvalue,3),") vs src mvalue",shape(mvalue)
        write(*,"(A,I2,A,I2)")       "Fortran DIMS: nssh_out(",nspecies,") vs src nssh",shape(nssh)
        write(*,"(A,I2,A,I2)")       "Fortran DIMS: nzx_out(",size(nzx,1),") vs src charges.nzx",shape(nzx)
        write(*,"(A,I2,A,I2)")       "Fortran DIMS: neighn_out(",natoms,") vs src neighn",shape(neighn)
        write(*,"(A,I2,A,I2,A,I2)")  "Fortran DIMS: neigh_j_out(",neigh_max,",",natoms,") vs src neigh_j",shape(neigh_j)
        write(*,"(A,I2,A,I2,A,I2)")  "Fortran DIMS: neigh_b_out(",neigh_max,",",natoms,") vs src neigh_b",shape(neigh_b)
        write(*,"(A,I2,A,I2,A,I2)")  "Fortran DIMS: xl_out(3,",size(xl,2),") vs src xl",shape(xl)
    end if
    h_mat_out = h_mat
    s_mat_out = s_mat

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

end subroutine firecore_get_HS_sparse

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

! Add other simple dimension getters if needed, e.g., for ME2c_max
subroutine firecore_get_ME2c_max(ME2c_max_out) bind(c, name='firecore_get_ME2c_max')
    use iso_c_binding
    use interactions, only: ME2c_max
    implicit none
    integer(c_int), intent(out) :: ME2c_max_out
    ME2c_max_out = ME2c_max
end subroutine firecore_get_ME2c_max
