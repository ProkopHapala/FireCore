

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
    use grid
    use debug
    use timing
    implicit none

    !interface
    !subroutine writeout_xsf (xsfname, message, aa)
    ! real, dimension (:), pointer, intent (in) :: aa
    ! character (len=40) xsfname
    ! character (len=30) message
    !end
    !end interface
    !interface
    !subroutine  project_orb(iband,ewfaux)
    !    integer iband
    !    real, dimension (:), pointer, intent (out) :: ewfaux
    !   end
    !end interface

    ! ====== global variables

    integer i,j, in1
    integer ikpoint
    integer imu
    real, dimension (3) :: k_temp
    real time_begin
    real time_end
    !real,  target,  dimension (:), allocatable :: ewfaux
    !real,  pointer, dimension (:)              :: pewf
    real,  dimension (:), allocatable :: ewfaux
    character (40) namewf
    character (30) mssg

    ! ====== Body

    iparam_file = 1
    !idebugWrite = 1
    idebugWrite = 0

    call cpu_time (time_begin)
    call initbasics ()    
    !call readdata ()
    call readdata_mcweda ()

    call init_wfs(norbitals, nkpoints)
    ikpoint = 1

    write (*,*) "fireball.f90: itheory, itheory_xc ", itheory, itheory_xc
    if(idebugWrite .gt. 0) write(*,*) "BEGIN fireball.f90 norbitals, nkpoints, max_scf_iterations", norbitals, nkpoints, max_scf_iterations
    ! TODO : It does not read CHARGES   (that is the reason for difference from Fireball-progs)
    write(*,*) "!!!! LOOP nstepf, max_scf_iterations ", nstepf, max_scf_iterations
    call init_FIRE( )

    call clean_ncall()
    do itime_step = 1,nstepf
        !debug_writeIntegral( interaction, isub, in1, in2, index )
        !max_scf_iterations = 1
        iforce = 0
        
        iforce = 1
        scf_achieved = .false.
        
        do Kscf = 1, max_scf_iterations
            !call clean_ncall()
            !write(*,*) "! ======== Kscf ", Kscf
            !call assemble_h ()
            call timer_start_i(1)
            call assemble_mcweda()
            call timer_stop_i(1)
            !call debug_writeBlockedMat( "S_mat.log", s_mat )
            !call debug_writeBlockedMat( "H_mat.log", h_mat )
            k_temp(:) = special_k(:,ikpoint)
            call timer_start_i(2)
            call solveH ( ikpoint, k_temp )
            call timer_stop_i(2)
            call timer_start_i(3)
            call denmat ()
            call timer_stop_i(3)
            sigma = sqrt(sum((Qin(:,:) - Qout(:,:))**2))
            write (*,*) "### SCF converged? ", scf_achieved, " Kscf ", Kscf, " |Qin-Oout| ",sigma," < tol ", sigmatol
            if( scf_achieved ) exit
            !Qin(:,:) = Qin(:,:)*(1.0-bmix) + Qout(:,:)*bmix   ! linear mixer 
            call mixer ()
            if(timing_verbosity .gt. 0) call writeclean_ncall()
        end do ! Kscf
        !call clean_ncall()
        call timer_start_i(4)
        !call postscf ()               ! optionally perform post-processing (DOS etc.)
        !call getenergy (itime_step)    ! calculate the total energy
        !call assemble_mcweda ()
        call getenergy_mcweda () 
        call timer_stop_i(4)
        call timer_start_i(11)
        call getforces_mcweda ()
        !call getforces ()   ! Assemble forces
        !do i=1, natoms
        !    write(*,*) "force[",i,"] ",  ftot(:,i)
        !end do
        call timer_stop_i(11)
        call move_ions_FIRE (itime_step, iwrtxyz )   ! Move ions now
        call write_bas()
        if(timing_verbosity .gt. 0) call writeclean_ncall()
        write (*,'(A,i6,A,i6,A,f16.8,A,f16.8)') " ###### Time Step", itime_step,"(of ",nstepf,") |Fmax| ",  deltaFmax , " > force_tol " , force_tol
        !	if ( FIRE_Ftot .lt. force_tol ) then
        if ( deltaFmax .lt. force_tol ) then
            write (*,*) ' +++++ FIRE.optionalimization converged +++++ '
            write (*,*) 'That`sall for now, bye ..'
            exit
        endif
    end do ! istepf

    if(iwrtewf .gt. 0) then
        !allocate   ( ewfaux(0:nrm-1))
        allocate   ( ewfaux(nrm))
        !pewf => ewfaux
        !call project_orb( 2, pewf )
        call project_orb( 2, ewfaux )
        write( namewf, '(A,i4.4,A)') 'bandplot_',1,'.xsf'
        mssg = 'density_3D'
        !call writeout_xsf ( namewf, mssg, pewf )
        call writeout_xsf ( namewf, mssg, ewfaux )
        deallocate ( ewfaux )
    end if

    call cpu_time (time_end)

    call writeclean_ncall()

    write (*,*) ' FIREBALL RUNTIME : ',time_end-time_begin,'[sec]'

    write(*,*) "TIME  1 assemble_mcweda  ", timer_sums(1)
    write(*,*) "TIME  2 solveH           ", timer_sums(2)
    write(*,*) "TIME  3 denmat           ", timer_sums(3)
    !write(*,*) "TIME  4 getenergy_mcweda ", timer_sums(4)
    write(*,*) "TIME 11 getforces_mcweda ", timer_sums(11)
    !write(*,*) "-----------------------------------"
    write(*,*) "TIME  5 assemble(Kscf=1)  ", timer_sums(5)
    !write(*,*) "TIME  6 get_ewald         ", timer_sums(6)
    !write(*,*) "TIME  7 assemble 1C       ", timer_sums(7)
    write(*,*) "TIME  8 assemble 2C       ", timer_sums(8)
    write(*,*) "TIME  9 assemble 3C       ", timer_sums(9)
    !write(*,*) "TIME 10 buildH            ", timer_sums(10)
    !write(*,*) "-----------------------------------"
    write(*,*) "TIME 20 assemle_2c(Kscf=1) ", timer_sums(20)
    write(*,*) "TIME 21 average_rho        ", timer_sums(21)
    !write(*,*) "TIME 22 assemble_olsxc_on  ", timer_sums(22)
    write(*,*) "TIME 23 assemble_olsxc_off ", timer_sums(23)
    write(*,*) "TIME 24 assemble_ca_2c     ", timer_sums(24)
    !write(*,*) "-----------------------------------"
    write(*,*) "TIME 31 assemble_3c        ", timer_sums(31)
    write(*,*) "TIME 32 assemble_3c_PP     ", timer_sums(32)
    write(*,*) "TIME 33 assemble_ca_3c     ", timer_sums(33)
    !write(*,*) "TIME 34 assemble_lr        ", timer_sums(34)
    !write(*,*) "-----------------------------------"
    !write(*,*) "TIME 50 assemble_F            ", timer_sums(50)
    write(*,*) "TIME 51 Dassemble_2c          ", timer_sums(51)
    !write(*,*) "TIME 52 Dassemble_2c_PP       ", timer_sums(52)
    write(*,*) "TIME 53 Dassemble_ca_olsxc_on ", timer_sums(53)
    write(*,*) "TIME 54 Dassemble_ca_olsxc_2c ", timer_sums(54)
    write(*,*) "TIME 55 Dassemble_ca_2c       ", timer_sums(55)
    write(*,*) "TIME 56 Dassemble_3c          ", timer_sums(56)
    write(*,*) "TIME 57 Dassemble_3c_PP       ", timer_sums(57)
    write(*,*) "TIME 58 Dassemble_ca_3c       ", timer_sums(58)
    write(*,*) "TIME 59 Dassemble_lr          ", timer_sums(59)
    write(*,*) "TIME 60 Dassemble_ca_olsxc_3c ", timer_sums(60)
    
    stop
end program fireball
