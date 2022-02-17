
! ==================================================
! ============ subroutine firecore_init
! ==================================================

! see : https://stackoverflow.com/questions/29153501/when-i-pass-a-pointer-from-fortran-to-c-the-value-of-the-first-element-is-lost
subroutine firecore_getPointer_ratom( ratom_ ) bind(c, name='firecore_getPointer_ratom')
    use iso_c_binding
    use configuration
    type(C_PTR),      intent(out) :: ratom_
    ratom_ = c_loc( ratom ) 
end subroutine

subroutine firecore_getPointer_ftot( ftot_ ) bind(c, name='firecore_getPointer_ftot')
    use iso_c_binding
    use forces
    type(C_PTR),      intent(out) :: ftot_
    ftot_ = c_loc( ftot ) 
end subroutine 

subroutine firecore_getPointer_charges( charges_ )  bind(c, name='firecore_getPointer_charges')
    use iso_c_binding
    use charges
    use options
    type(C_PTR),      intent(out) :: charges_
    if (iqout .eq. 2) then
        charges_ = c_loc( QMulliken_TOT )
    else
        charges_ = c_loc( QLowdin_TOT   )
    end if
    !write(*,*) "DEBUG firecore_getCharges END"
end 

subroutine firecore_init( natoms_, atomTypes, atomsPos ) bind(c, name='firecore_init')
    use iso_c_binding
    use options
    use loops
    !use fire
    use configuration
    !use energy
    !use forces
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
    ! ====== global variables
    !real time_begin
    !real time_end
    integer i, ispec, in1, iatom, numorbPP_max
    real distance
    real, dimension (3) :: vector
    logical zindata
    ! ====== Body
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
    write(*,*) "DEBUG initbasics START "
    call initconstants ! (sigma, sigmaold, scf_achieved)
    write(*,*) "DEBUG 1 natoms, natoms_", natoms, natoms_
    scf_achieved = .true.
    write(*,*) "DEBUG 2 natoms, natoms_", natoms, natoms_
    call diagnostics (ioff2c, ioff3c, itestrange, testrange)    ! IF_DEF_DIAGNOSTICS
    write(*,*) "DEBUG 3 natoms, natoms_", natoms, natoms_
    call readparam ()
    write(*,*) "DEBUG 4 natoms, natoms_", natoms, natoms_

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
    write(*,*) "DEBUG shape(iatyp) ",shape(iatyp),"DEBUG shape(atomTypes) ", shape(atomTypes)
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
    call write_to_xyz( "#DEBUG libFireCore::firecore_init() ", 1 )
    return
end subroutine firecore_init

! ==================================================
! ============ subroutine firecore_SCF
! ==================================================

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
end 

! ==================================================
! ============ subroutine init_wfs
! ==================================================

