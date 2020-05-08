

program fireball
    use options
    use loops
    use interactions
    use integrals
    use density
    use kpoints
    use charges
    
    implicit none

    ! ====== global variables

    integer i,j
    integer ikpoint
    integer imu
    real, dimension (3) :: k_temp

    real time_begin
    real time_end

    ! ====== Body
    call cpu_time (time_begin)
    write(*,*) "DEBUG 1"
    call initbasics ()

    do  i = 1,23
        write (*,'(A,i10)',advance='no') "DEBUG ind2c ", i
        do j = 0,8
            write (*,'(i10)',advance='no')   ind2c( i,j)
        end do
        write (*,*)
    end do

    ! TODO : It does not read CHARGES   (that is the reason for difference from Fireball-progs)

    write(*,*) "DEBUG 2"
    !call readdata ()
    call readdata_mcweda ()
    write(*,*) "DEBUG 3"
    !call main_loop ()
    !call scf_loop (itime_step)

    ! =========== Allocate MOs
    nkpoints = 1
    if (iqout .ne. 2) then 
        allocate (blowre (norbitals, norbitals, nkpoints))
        allocate (blowim (norbitals, norbitals, nkpoints))
        blowre = 0.0d0
        blowim = 0.0d0
    endif
    allocate (bbnkre (norbitals, norbitals, nkpoints))
    allocate (bbnkim (norbitals, norbitals, nkpoints))
    allocate (eigen_k (norbitals, nkpoints))
    eigen_k = 0.0d0
    bbnkim = 0.0d0
    bbnkre = 0.0d0
    write(*,*) "DEBUG norbitals, nkpoints ", norbitals, nkpoints


    do Kscf = 1, max_scf_iterations
        write(*,*) "! ======== Kscf ", Kscf
        call assemble_h ()
        !call diag_k()
        k_temp(:) = special_k(:,ikpoint)
        !call kspace( nprocs, my_proc, Kscf, iqout, icluster, iwrteigen, ikpoint, k_temp, nkpoints, iwrtdos, iwrthop, iwrtatom, itrans, igap )
        call solveH   ( Kscf, ikpoint, k_temp )
        !call build_rho() 
        call denmat ()
        Qin(:,:) = Qin(:,:)*(1.0-bmix) + Qout(:,:)*bmix   ! linear mixer 
        !call mixCharge
    end do
    write(*,*) "DEBUG 4"
    !call postscf ()               ! optionally perform post-processing (DOS etc.)
    !call getenergy (itime_step)    ! calculate the total energy
    call getenergy_mcweda () 
    write(*,*) "DEBUG 5"

    call cpu_time (time_end)
    write (*,*) ' FIREBALL RUNTIME : ',time_end-time_begin,'[sec]'

    stop
end program fireball
