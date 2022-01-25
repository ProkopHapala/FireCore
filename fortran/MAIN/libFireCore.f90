
! ==================================================
! ============ subroutine firecore_init
! ==================================================

subroutine firecore_init( natoms_, atomTypes, atomsPos )
    use options, only : idebugWrite
    use loops
    !use fire
    use configuration
    !use energy
    !use forces
    use interactions
    !use integrals
    !use density
    use kpoints
    !use charges
    !use debug
    implicit none
    ! ====== Parameters
    integer,                      intent(in) :: natoms_
    integer, dimension(:),        intent(in) :: atomTypes
    real,    dimension(3,natoms), intent(in) :: atomsPos
    ! ====== global variables
    !real time_begin
    !real time_end
    ! ====== Body
    idebugWrite = 0
    !call cpu_time (time_begin)
    call initbasics ()    
    !call readdata ()
    call readdata_mcweda ()
    call init_wfs(norbitals, nkpoints)
    !if(idebugWrite .gt. 0) write(*,*) "DEBUG fireball.f90 norbitals, nkpoints, max_scf_iterations", norbitals, nkpoints
    ! TODO : It does not read CHARGES   (that is the reason for difference from Fireball-progs)
    !write(*,*) "!!!! LOOP nstepf, max_scf_iterations ", nstepf, max_scf_iterations
    !call init_FIRE( )
    !call cpu_time (time_end)
    !write (*,*) ' FIREBALL RUNTIME : ',time_end-time_begin,'[sec]'
    stop
end subroutine firecore_init

! ==================================================
! ============ subroutine firecore_SCF
! ==================================================

subroutine firecore_SCF( nmax_scf, forces_ )
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
    integer,                 intent(in) :: nmax_scf
    real,    dimension(:,:), intent(in) :: forces_
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
    !call cpu_time (time_end)
    !write (*,*) ' FIREBALL RUNTIME : ',time_end-time_begin,'[sec]'
    stop
end subroutine firecore_SCF

! ==================================================
! ============ subroutine init_wfs
! ==================================================

