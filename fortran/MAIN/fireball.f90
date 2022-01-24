

program fireball
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

    ! ====== global variables

    integer i,j, in1
    integer ikpoint
    integer imu
    real, dimension (3) :: k_temp

    real time_begin
    real time_end

    ! ====== Body

    !idebugWrite = 1
    idebugWrite = 0

    call cpu_time (time_begin)
    call initbasics ()    
    !call readdata ()
    call readdata_mcweda ()

    call init_wfs(norbitals, nkpoints)
    ikpoint = 1

    if(idebugWrite .gt. 0) write(*,*) "DEBUG fireball.f90 norbitals, nkpoints, max_scf_iterations", norbitals, nkpoints, max_scf_iterations
    ! TODO : It does not read CHARGES   (that is the reason for difference from Fireball-progs)
    write(*,*) "!!!! LOOP nstepf, max_scf_iterations ", nstepf, max_scf_iterations
    call init_FIRE( )

    do itime_step = 1,nstepf
        !debug_writeIntegral( interaction, isub, in1, in2, index )
        !max_scf_iterations = 1
        iforce = 0
        
        iforce = 1
        scf_achieved = .false.
        
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
        !do i=1, natoms
        !    write(*,*) "force[",i,"] ",  ftot(:,i)
        !end do

        call move_ions_FIRE (itime_step, iwrtxyz )   ! Move ions now
        call write_bas()

        write (*,'(A,i6,A,i6,A,f16.8,A,f16.8)') " ###### Time Step", itime_step,"(of ",nstepf,") |Fmax| ",  deltaFmax , " > force_tol " , force_tol
        !	if ( FIRE_Ftot .lt. force_tol ) then
        if ( deltaFmax .lt. force_tol ) then
            write (*,*) ' +++++ FIRE.optionalimization converged +++++ '
            write (*,*) 'That`sall for now, bye ..'
            exit
        endif
    end do ! istepf

    call cpu_time (time_end)
    write (*,*) ' FIREBALL RUNTIME : ',time_end-time_begin,'[sec]'

    stop
end program fireball

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
    write(*,*) "init_wfs norb, nkpt ",  norb, nkpt, nkpoints
    allocate (blowre (norb, norb, nkpoints))
    allocate (blowim (norb, norb, nkpoints))
    blowre = 0.0d0
    blowim = 0.0d0
    !endif
    allocate (bbnkre (norb, norb, nkpoints))
    allocate (bbnkim (norb, norb, nkpoints))
    allocate (eigen_k (norb, nkpoints))
    write(*,*) "shape(eigen_k)", shape(eigen_k) 
    eigen_k = 0.0d0
    bbnkim = 0.0d0
    bbnkre = 0.0d0

    ! --- Kpoints
    allocate (special_k(3,nkpoints))
    allocate (weight_k (nkpoints))
    special_k(:,:) = 0
    weight_k (:)   = 1

end subroutine init_wfs