

program fireball
    use options
    use loops
    use configuration
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

    idebugWrite = 1
    !idebugWrite = 0

    call cpu_time (time_begin)
    call initbasics ()

    ! TODO : It does not read CHARGES   (that is the reason for difference from Fireball-progs)

    !call readdata ()
    call readdata_mcweda ()
    !call main_loop ()
    !call scf_loop (itime_step)

    ! =========== Allocate MOs
    nkpoints = 1
    ! k_temp(:) = 0
    if (iqout .ne. 2) then 
        allocate (blowre (norbitals, norbitals, nkpoints))
        allocate (blowim (norbitals, norbitals, nkpoints))
        blowre = 0.0d0
        blowim = 0.0d0
    endif
    allocate (bbnkre (norbitals, norbitals, nkpoints))
    allocate (bbnkim (norbitals, norbitals, nkpoints))
    allocate (eigen_k (norbitals, nkpoints))
    allocate (special_k(3,nkpoints))
    allocate (weight_k (nkpoints))
    eigen_k = 0.0d0
    bbnkim = 0.0d0
    bbnkre = 0.0d0
    if(idebugWrite .gt. 0) write(*,*) "DEBUG fireball.f90 norbitals, nkpoints, max_scf_iterations", norbitals, nkpoints, max_scf_iterations

    ikpoint = 1
    special_k(:,:) = 0
    weight_k (:)   = 1
    k_temp(:) = special_k(:,ikpoint)


    !debug_writeIntegral( interaction, isub, in1, in2, index )
    max_scf_iterations = 10
    iforce = 0
    do Kscf = 1, max_scf_iterations
        write(*,*) "! ======== Kscf ", Kscf
        !call assemble_h ()

        call assemble_mcweda ()
        !call diag_k()
        !write(*,*) "shape(k_temp) ", shape(k_temp)
        !write(*,*) "special_k     ", shape(special_k)
        !k_temp(:) = special_k(:,ikpoint)
        !call kspace( nprocs, my_proc, Kscf, iqout, icluster, iwrteigen, ikpoint, k_temp, nkpoints, iwrtdos, iwrthop, iwrtatom, itrans, igap )

        !call debug_writeBlockedMat( "S_mat.log", s_mat )
        !call debug_writeBlockedMat( "H_mat.log", h_mat )
        !write (*,*) "DEBUG STOP in fireball.f90"
        !stop  ! DEBUG

        call solveH   ( ikpoint, k_temp )
        !call build_rho() 
        call denmat ()
        
        !write (*,*) "Qin ",  Qin(1,:)
        !write (*,*) "Qout ", Qout(1,:)

        sigma = sqrt(sum((Qin(:,:) - Qout(:,:))**2))

        if ( sigma .lt. sigmatol) then
            write (*,*) "# SCF converged ", Kscf ,sigma, sigmatol
            exit
        else 
            write (*,*) "# SCF converged not ", Kscf ,sigma, sigmatol
        end if ! simga

        !Qin(:,:) = Qin(:,:)*(1.0-bmix) + Qout(:,:)*bmix   ! linear mixer 
        !call mixCharge
        call mixer ()

        write (*,*) "Qin ",  Qin(1,:)
        write (*,*) "Qout ", Qout(1,:)

        !do i = 1, natoms
        !    in1 = imass(i)
        !    write (*,'(A,i5,A)',advance='no') "atom[",i,"]" 
        !    do j = 1, nssh(in1)
        !      write (*,'(f10.5)',advance='no') Qin (j,i)
        !      write (*,'(f10.5)',advance='no') Qout(j,i)
        !    end do
        !    write (*,*)
        !end do ! 

    end do ! Kscf
    !call postscf ()               ! optionally perform post-processing (DOS etc.)
    !call getenergy (itime_step)    ! calculate the total energy
    iforce = 1
    call assemble_mcweda ()
    call getenergy_mcweda () 
    call getforces ()   ! Assemble forces
    
    !call move_ions_FIRE (itime_step, 0 )   ! Move ions now

    do i=1, natoms
        write(*,*) "force[",i,"] ",  ftot(:,i)
    end do

    call cpu_time (time_end)
    write (*,*) ' FIREBALL RUNTIME : ',time_end-time_begin,'[sec]'

    stop
end program fireball
