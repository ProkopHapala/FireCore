
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
    if (iqout .eq. 2) then
        charges_ = c_loc( QMulliken_TOT )
    else
        charges_ = c_loc( QLowdin_TOT   )
    end if
    !write(*,*) "DEBUG firecore_getCharges END"
end subroutine

subroutine firecore_preinit( )  bind(c, name='firecore_preinit' )
    use iso_c_binding
    use options
    use configuration
    use kpoints
    !use MD
    use loops
    use charges
    !use barrier
    !use nonadiabatic
    use integrals !, only : fdataLocation
    implicit none
    ! ========= body
    iparam_file = 0
    call initconstants ! (sigma, sigmaold, scf_achieved)
    scf_achieved = .true.
    call diagnostics (ioff2c, ioff3c, itestrange, testrange)    ! IF_DEF_DIAGNOSTICS
    !call readparam ()
    call set_default_params ()
    igrid = 1
    call checksum_options ()
    write(*,*) "DEBUG firecore_preinit() icluster", icluster
end subroutine

subroutine firecore_set_lvs( lvs_ )  bind(c, name='firecore_set_lvs')
    use iso_c_binding
    use configuration
    use options
    implicit none
    real(c_double), dimension(3,3), intent(in) :: lvs_
    icluster = 0
    a1vec = lvs_(:,1)
    a2vec = lvs_(:,2)
    a3vec = lvs_(:,3)
    write(*,*) "DEBUG firecore_set_lvs() icluster", icluster
end subroutine

subroutine firecore_init( natoms_, atomTypes, atomsPos ) bind(c, name='firecore_init')
    use iso_c_binding
    use options
    use loops
    !use fire
    use configuration
    !use energy
    !use forces
    use grid
    use interactions
    use integrals
    use density
    use kpoints
    use charges
    use FIRE, only: write_to_xyz
    !use debug
    implicit none

    ! ====== Parameters
    integer(c_int),                      intent(in), value :: natoms_
    !integer(c_int), dimension(natoms),   intent(in)        :: atomTypes
    !real(c_double), dimension(3,natoms), intent(in)        :: atomsPos
    integer(c_int), dimension(natoms_),   intent(in)        :: atomTypes
    real(c_double), dimension(3,natoms_), intent(in)        :: atomsPos
    !integer(c_int),                      intent(out)        :: norb
    ! ====== global variables
    !real time_begin
    !real time_end
    integer i, ispec, in1, iatom, numorbPP_max
    real distance
    real, dimension (3) :: vector
    logical zindata
    ! ====== Body
    write(*,*) "DEBUG firecore_init() icluster", icluster
    natoms = natoms_

    write(*,*) "DEBUG natoms, natoms_", natoms, natoms_
    write(*,*) "DEBUG shape(atomTypes,atomsPos)", shape(atomTypes), shape(atomsPos)
    do i = 1,natoms_
        write(*,*) atomTypes(i), atomsPos(:,i)
    end do !

    idebugWrite = 0
    !call cpu_time (time_begin)

    !call initbasics () 
    ! >>> BEGIN INIT_BASICS
    !call initconstants ! (sigma, sigmaold, scf_achieved)
    !scf_achieved = .true.
    !call diagnostics (ioff2c, ioff3c, itestrange, testrange)    ! IF_DEF_DIAGNOSTICS
    !call readparam ()

    write(*,*) "icluster", icluster
    write(*,*) "a1vec", a1vec
    write(*,*) "a3vec", a2vec
    write(*,*) "a2vec", a3vec
    write(*,*) "a2vec", a3vec
    write(*,*) "rm1,rm2,rm3,nrm ", rm1,rm2,rm3, nrm
    write(*,*) "em1,em2,em3,nem ", rm1,rm2,rm3, nem


 ! Allocate more arrays.
    allocate (degelec (natoms))
    allocate (iatyp (natoms))
    allocate (imass (natoms))
    allocate (ratom (3, natoms))
    allocate (nowMinusInitialPos (3, natoms))
    allocate (initialPosition (3, natoms))
    allocate (vatom (3, natoms))
    allocate (symbol (natoms))
    allocate (xmass (natoms))
    allocate (ximage (3, natoms))
    ! allocate (mask (3,natoms))      ! IF_DEF_OPT_END
    ! mask = 1.0d0                    ! IF_DEF_OPT_END
    ximage = 0.0d0
    !call readbasis (nzx, imass)
    !write(*,*) "DEBUG shape(iatyp) ",shape(iatyp),"DEBUG shape(atomTypes) ", shape(atomTypes)
    iatyp(:) = atomTypes(:)
    ratom(:,:) = atomsPos(:,:)

 ! Read the info.dat file.  Allocate lsshPP, etc. inside!
    call readinfo ()

    write(*,*) "DEBUG natoms,natoms_ ", natoms,natoms_
    write(*,*) "DEBUG shape(ratom) ",shape(ratom),"DEBUG shape(atomsPos) ", shape(atomsPos)
    !imass(:)   = atomTypes(:)
    
    call cross (a2vec, a3vec, vector)
    Vouc=a1vec(1)*vector(1)+a1vec(2)*vector(2)+a1vec(3)*vector(3)
    call initboxes (1)

    !write(*,*) "DEBUG nssh(:)", nssh(:)

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
    write (*,*) ' Initiallizing arrays '
    call allocate_neigh()
    call allocate_f() 
    call allocate_h() 
    call allocate_rho() 
    !call allocate_dos() ! IF_DEF_GRID_DOS
    !<<< END INIT_BASICS
    !write(*,*) "DEBUG initbasics END "

    !call readdata ()
    call readdata_mcweda ()
    !write(*,*) "DEBUG readdata_mcweda END "
    call init_wfs(norbitals, nkpoints)
    !write(*,*) "DEBUG firecore_init END "
    write (*,*) ", norbitals ", norbitals
    !norb = norbitals
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
    Kscf       = Kscf_
    iforce     = iforce_
    ratom(:,:) = positions_(:,:)
    idebugWrite = 0
    verbosity   = 0 
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
    ratom(:,:) = positions_(:,:)
    !idebugWrite = 0
    !verbosity   = 0 
    iforce      = iforce_
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



subroutine firecore_evalForce( nmax_scf, positions_, forces_ )  bind(c, name='firecore_evalForce')
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
    !real(c_double), dimension(:,:), intent(out)      :: forces_
    real(c_double), dimension(3,natoms), intent(in) :: positions_
    real(c_double), dimension(3,natoms), intent(out) :: forces_
    ! ====== global variables
    integer ikpoint, i
    real, dimension (3) :: k_temp
    !real time_begin
    !real time_end
    ! ====== Body

    ratom(:,:) = positions_(:,:)
    !idebugWrite = 1
    !idebugWrite = 0
    !call cpu_time (time_begin)

    idebugWrite = 0
    verbosity   = 0 
    iforce      = 1

    ftot(:,:) = 0
    ikpoint = 1
    scf_achieved = .false.
    max_scf_iterations = nmax_scf
    if(verbosity.gt.0)write(*,*) "!!!! SCF LOOP max_scf_iterations ", max_scf_iterations, scf_achieved
    do Kscf = 1, max_scf_iterations
        if(idebugWrite.gt.0)write(*,*) "! ======== Kscf ", Kscf
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
    !call postscf ()               ! optionally perform post-processing (DOS etc.)
    !call getenergy (itime_step)    ! calculate the total energy
    !call assemble_mcweda ()
    call getenergy_mcweda ()
    call getforces_mcweda ()
    !do i = 1, natoms
    !    write (*,*) "Force[",i,"] ", ftot(:,i)
    !end do
    !call getforces ()   ! Assemble forces
    forces_(:,:) = ftot(:,:)
    !call cpu_time (time_end)
    !write (*,*) ' FIREBALL RUNTIME : ',time_end-time_begin,'[sec]'
    if(verbosity.gt.0)write (*,*) '!!!! SCF LOOP DONE in ', Kscf, " iterations"
    return
end subroutine firecore_evalForce

subroutine firecore_getCharges( charges_ )  bind(c, name='firecore_getCharges')
    use iso_c_binding
    use configuration
    use charges
    use options
    real(c_double), dimension(natoms), intent(out) :: charges_
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
    write (*,*) "shape(bbnkre) ", shape(bbnkre), ikp
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
    write (*,*) "shape(bbnkre) ", shape(bbnkre), ikp, iMO
    bbnkre(:,iMO,ikp) = wfcoefs(:)
end subroutine

subroutine firecore_setupGrid( Ecut_, ifixg0_, g0_,    ngrid, dCell  )  bind(c, name='firecore_setupGrid' )
    use iso_c_binding
    use grid
    use configuration
    implicit none
    ! ========= Parameters
    !integer(c_int),                 intent(in),value :: nmax_scf
    !real(c_double), dimension(3,natoms), intent(in) :: positions_
    !real(c_double), dimension(3,natoms), intent(out) :: forces_
    real   (c_double)               ,intent(in),value :: Ecut_
    integer(c_int)                  ,intent(in),value :: ifixg0_
    real   (c_double), dimension (3),intent(in) :: g0_
    integer(c_int),    dimension (3),intent  (out) :: ngrid
    real   (c_double), dimension (3,3),intent(out) :: dCell
    !call readgrid !(iwrtewf)
    ! Namelist /mesh/ Ecut, iewform, npbands, pbands, ewfewin_max, ewfewin_min, ifixg0, g0
    ! ========= Body
    Ecut   = Ecut_
    ifixg0 = ifixg0_
    g0(:)  = g0_(:)
    call allocate_grid !(natoms, nspecies)
    call read_wf ()
    call read_vna ()
    ! np(i) = int (2 * sqrt(Ecut) / (cvec(i)*abohr) + 1)     
    call initgrid !(icluster)
    ngrid(1)=rm1
    ngrid(2)=rm2
    ngrid(3)=rm3
    dCell(:,:) = elvec(:,:)
end subroutine

subroutine firecore_getGridMO( iMO, ewfaux )  bind(c, name='firecore_getGridMO' )
    use iso_c_binding
    use grid
    use configuration
    implicit none
    ! ========= Parameters
    integer(c_int), value :: iMO
    real(c_double), dimension (nrm), intent(out) :: ewfaux
    ! ========= Body
    !allocate   ( ewfaux(0:nrm-1))
    !pewf => ewfaux
    !write(*,*) "firecore_getGridMO ", iMO
    call project_orb( iMO, ewfaux )
end subroutine

subroutine firecore_getGridDens( ewfaux )  bind(c, name='firecore_getGridDens' )
    use iso_c_binding
    use grid
    use configuration
    implicit none
    ! ========= Parameters
    !integer(c_int), value :: imo0, imo1
    real(c_double), dimension (nrm), intent(out) :: ewfaux
    ! ========= Body
    !write(*,*) "firecore_getGridDens "
    call project_dens( ewfaux )
end subroutine

subroutine firecore_orb2xsf( iMO )  bind(c, name='firecore_orb2xsf' )
    use iso_c_binding
    use options, only: icluster
    use grid
    implicit none
 ! Local Parameters and Data Declaration
 ! ===========================================================================
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
 ! Local Parameters and Data Declaration
    real(c_double), value :: f_den0
    ! ========= Variables
    real, dimension (:), allocatable :: ewfaux
    character(40)   :: namewf
    character(4)    :: name
    character (len=30) mssg
    integer i
    ! ========= Body
    allocate ( ewfaux(0:nrm-1))
    call project_dens( ewfaux )
    if( f_den0*f_den0 > 1.e-16 ) then
        call project_dens0( f_den0, ewfaux )
    end if 
    namewf = 'density.xsf'
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
    do i=1,n
        call getpsi(in1,issh,NORM2(poss(:,i)),psi,dpsi)
        call getYlm(l,poss(:,i),Y,dY) 
        ys(i) = psi  *Y(m)
        !write (*,*) "", in1,issh,l,m, i,x,psi,Y(m)
    end do 
end subroutine

! ==================================================
! ============ subroutine init_wfs
! ==================================================

