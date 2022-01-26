
! ==================================================
! ============ subroutine firecore_init
! ==================================================

subroutine firecore_hello( ) bind(c, name='firecore_hello')
    use iso_c_binding
    ! ====== Parameters
    ! ====== global variables
    
    ! ====== Body
    write(*,*) " firecore_hello says HELLO !!!! "
    return
end subroutine firecore_hello

function sum2(a) result(b) bind(c, name='sum2')
    use iso_c_binding
    implicit none
    real(c_double), intent(in)  :: a
    real(c_double)              :: b
    b = a + 2.d0
end function sum2

function sum2val(a,b) result(c) bind(c, name='sum2val')
    use iso_c_binding
    implicit none
    real(c_double), intent(in),value  :: a
    real(c_double), intent(in),value  :: b
    real(c_double)                    :: c
    c = a + b
end function sum2val

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
    integer i, ispec, iatom, numorbPP_max
    real distance
    real, dimension (3) :: vector
    logical zindata
    ! ====== Body
    natoms = natoms_

    write(*,*) "DEBUG natoms_", natoms
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
    scf_achieved = .true.
    call diagnostics (ioff2c, ioff3c, itestrange, testrange)    ! IF_DEF_DIAGNOSTICS
    call readparam ()
    call readinfo () 
    allocate (degelec   (natoms))
    allocate (imass     (natoms))
    allocate (symbol    (natoms))
    allocate (xmass     (natoms))
    allocate (ratom  (3, natoms))
    allocate (vatom  (3, natoms))
    allocate (ximage (3, natoms))
    allocate (nowMinusInitialPos (3, natoms))
    allocate (initialPosition    (3, natoms))
    ximage = 0.0d0
    !call readbasis (nzx, imass)

    !imass(:)   = atomTypes(:)
    ratom(:,:) = atomsPos (:,:)
    write (*,*) " check atom species ... "
    do iatom = 1, natoms
        zindata = .false.
        do ispec = 1, nspecies
             if ( atomTypes(iatom) .eq. nzx(ispec)) then
                zindata = .true.
                imass(:) = ispec
             end if
        end do
        if ( .not. zindata ) then
            write (*,*) "atom[",iatom,"] type not know ", atomTypes(iatom)
            stop
            return
        else 
            write (*,*) "atom[",iatom,"] ",imass(iatom) 
        end if ! zindata

    end do ! iatom
    write(*,*) " ... Atoms are OK "
    !return

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
    write (*,*) ' Initiallizing arrays '
    call allocate_neigh()
    call allocate_f() 
    call allocate_h() 
    call allocate_rho() 
    !call allocate_dos() ! IF_DEF_GRID_DOS
    !<<< END INIT_BASICS
    write(*,*) "DEBUG initbasics END "

    !call readdata ()
    call readdata_mcweda ()
    write(*,*) "DEBUG readdata_mcweda END "
    call init_wfs(norbitals, nkpoints)
    write(*,*) "DEBUG firecore_init END "
    return
end subroutine firecore_init

! ==================================================
! ============ subroutine firecore_SCF
! ==================================================

subroutine firecore_evalForce( nmax_scf, forces_ )  bind(c, name='firecore_evalForce')
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
    real(c_double), dimension(:,:), intent(out)      :: forces_
    ! ====== global variables
    integer ikpoint
    real, dimension (3) :: k_temp
    !real time_begin
    !real time_end
    ! ====== Body
    !idebugWrite = 1
    !idebugWrite = 0
    !call cpu_time (time_begin)
    ikpoint = 1
    max_scf_iterations = nmax_scf
    write(*,*) "!!!! SCF LOOP max_scf_iterations ", max_scf_iterations
    do Kscf = 1, max_scf_iterations
        !write(*,*) "! ======== Kscf ", Kscf
        !call assemble_h ()
        call assemble_mcweda ()
        !call debug_writeBlockedMat( "S_mat.log", s_mat )
        !call debug_writeBlockedMat( "H_mat.log", h_mat )
        k_temp(:) = special_k(:,ikpoint)
        call solveH ( ikpoint, k_temp )
        call denmat ()
        sigma = sqrt(sum((Qin(:,:) - Qout(:,:))**2))
        write (*,*) "### SCF converged? ", scf_achieved, " Kscf ", Kscf, " |Qin-Oout| ",sigma," < tol ", sigmatol
        if( scf_achieved ) exit
        !Qin(:,:) = Qin(:,:)*(1.0-bmix) + Qout(:,:)*bmix   ! linear mixer 
        call mixer ()
    end do ! Kscf
    !call postscf ()               ! optionally perform post-processing (DOS etc.)
    !call getenergy (itime_step)    ! calculate the total energy
    !call assemble_mcweda ()
    call getenergy_mcweda () 
    call getforces ()   ! Assemble forces

    forces_(:,:) = ftot(:,:)
    !call cpu_time (time_end)
    !write (*,*) ' FIREBALL RUNTIME : ',time_end-time_begin,'[sec]'
    return
end subroutine firecore_evalForce

! ==================================================
! ============ subroutine init_wfs
! ==================================================

