
module timing
    implicit none

    integer, parameter :: timer_index_max = 100
    integer, parameter :: max_interaction = 23
    integer timer_index

    integer ncall_doscentros
    integer ncall_doscentrosS
    integer ncall_trescentros
    integer ncall_trescentrosS
    integer ncall_Dtrescentros
    integer ncall_DtrescentrosS

    integer ncall_epsilon
    integer ncall_deps2cent
    integer ncall_interpolate_1d
    integer ncall_interpolate_2d
    integer ncall_twister
    integer ncall_twisterd
    integer ncall_makeDmat
    integer ncall_rotate_fb
    integer ncall_interpolate_2d_vec
    integer ncall_interpolate_1d_vec



    integer ncall_assemble_sVNL
    integer ncall_assemble_2c
    integer ncall_assemble_2c_PP
    integer ncall_average_ca_rho
    integer ncall_average_rho
    integer ncall_assemble_snxc_on
    integer ncall_assemble_snxc_off
    integer ncall_assemble_olsxc_on
    integer ncall_assemble_olsxc_off
    integer ncall_assemble_ca_2c
    integer ncall_assemble_3c
    integer ncall_assemble_3c_PP
    integer ncall_assemble_ca_3c
    integer ncall_assemble_lr

!           interaction    subtype
!               1            isorp = 0, 8    average density 
!               1            0               overlap
!               2            0               vna_ontopl
!               2            isorp = 1, 8    vna_ontopl shell isorp
!               3            0               vna_ontopr
!               3            isorp = 1, 8    vna_ontopr shell isorp
!               4            0               vna_atom
!               4            isorp = 1, 8    vna_atom  shell isorp
!               5            0               non-local
!               6            0, 4            xc_ontop
!               7            0, 4            xc_atom
!               8            0, 4            xc_correction
!               9            0               z-dipole
!               10           0               y-dipole
!               11           0               x-dipole
!               12           0               coulomb
!               13           0               kinetic
!               14           0               extended-hubbard
!               15           isorp = 1, 8    density_ontopl
!               16           isorp = 1, 8    density_ontopr
!               17           isorp = 1, 8    density_atom       !  spherical density (only doscentrosS)
!               18           isorp = 1, 4    sph dens ontopl
!               19           isorp = 1, 4    sph dens_ontopr
!               20           isorp = 1, 4    sph den_atom
!               21           isorp = 0       sph overlap

!    interaction .eq. 1  bcna_01
!    interaction .eq. 2  xc3c_01
!    interaction .eq. 3  den3_01
!    interaction .eq. 4  den3S_01 

   integer ncall_interpolate_1d_inter(max_interaction)
   integer ncall_interpolate_2d_inter(max_interaction)

    

    ! integer ntotcall_doscentros
    ! integer ntotcall_doscentrosS
    ! integer ntotcall_trescentros
    ! integer ntotcall_trescentrosS
    ! integer ntotcall_epsilon
    ! integer ntotcall_deps2cent
    ! integer ntotcall_interpolate_1d
    ! integer ntotcall_interpolate_2d
    ! integer ntotcall_twister
    ! integer ntotcall_twisterd
    ! integer ntotcall_makeDmat
    ! integer ntotcall_rotate_fb
    ! integer ntotcall_interpolate_2d_vec

    real timer_starts(timer_index_max)
    real timer_sums  (timer_index_max)

    integer interaction_glob

    contains

    ! subroutine clean_ntotcall()
    !     ntotcall_doscentros   =0
    !     ntotcall_doscentrosS  =0
    !     ntotcall_trescentros  =0
    !     ntotcall_trescentrosS =0
    !     ntotcall_epsilon      =0
    !     ntotcall_deps2cent    =0
    !     ntotcall_interpolate_1d =0
    !     ntotcall_interpolate_2d =0
    !     ntotcall_twister      =0
    !     ntotcall_twisterd     =0
    !     ntotcall_makeDmat     =0
    !     ntotcall_rotate_fb    =0
    !     ntotcall_interpolate_2d_vec = 0
    ! end subroutine clean_ncall

    subroutine clean_ncall()
        ncall_doscentros    =0
        ncall_doscentrosS   =0
        ncall_trescentros   =0
        ncall_trescentrosS  =0
        ncall_Dtrescentros  =0
        ncall_DtrescentrosS =0
        ncall_epsilon       =0
        ncall_deps2cent     =0
        ncall_interpolate_1d=0
        ncall_interpolate_2d=0
        ncall_twister       =0
        ncall_twisterd      =0
        ncall_makeDmat      =0
        ncall_rotate_fb     =0
        ncall_interpolate_2d_vec = 0
        ncall_interpolate_1d_vec = 0

        ncall_interpolate_1d_inter(:) = 0
        ncall_interpolate_2d_inter(:) = 0

        ncall_assemble_sVNL      = 0
        ncall_assemble_2c        = 0
        ncall_assemble_2c_PP     = 0
        ncall_average_ca_rho     = 0
        ncall_average_rho        = 0
        ncall_assemble_snxc_on   = 0
        ncall_assemble_snxc_off  = 0
        ncall_assemble_olsxc_on  = 0
        ncall_assemble_olsxc_off = 0
        ncall_assemble_ca_2c     = 0
        ncall_assemble_3c        = 0
        ncall_assemble_3c_PP     = 0
        ncall_assemble_ca_3c     = 0
        ncall_assemble_lr        = 0

        ! ntotcall_doscentros   = ntotcall_doscentros 
        ! ntotcall_doscentrosS  = ntotcall_doscentrosS
        ! ntotcall_trescentros  = ntotcall_trescentros
        ! ntotcall_trescentrosS = ntotcall_trescentrosS
        ! ntotcall_epsilon      = ntotcall_epsilon
        ! ntotcall_deps2cent    = ntotcall_deps2cent
        ! ntotcall_interpolate_1d = ntotcall_interpolate_1d
        ! ntotcall_interpolate_2d = ntotcall_interpolate_2d
        ! ntotcall_twister      = ntotcall_twister 
        ! ntotcall_twisterd     = ntotcall_twisterd 
        ! ntotcall_makeDmat     = ntotcall_makeDmat
        ! ntotcall_rotate_fb    = ntotcall_rotate_fb 
        ! ntotcall_interpolate_2d_vec = ntotcall_interpolate_2d_vec
    end subroutine clean_ncall

    subroutine write_ncall()
        integer i
        write(*,*) "ncall_doscentros         ", ncall_doscentros
        write(*,*) "ncall_doscentrosS        ", ncall_doscentrosS  
        write(*,*) "ncall_trescentros        ", ncall_trescentros  
        write(*,*) "ncall_trescentrosS       ", ncall_trescentrosS 
        write(*,*) "ncall_Dtrescentros       ", ncall_Dtrescentros  
        write(*,*) "ncall_DtrescentrosS      ", ncall_DtrescentrosS 
        write(*,*) "ncall_epsilon            ", ncall_epsilon      
        write(*,*) "ncall_deps2cent          ", ncall_deps2cent    
        write(*,*) "ncall_interpolate_1d     ", ncall_interpolate_1d 
        write(*,*) "ncall_interpolate_1d_vec ", ncall_interpolate_1d_vec
        write(*,*) "ncall_interpolate_2d     ", ncall_interpolate_2d 
        write(*,*) "ncall_interpolate_2d_vec ", ncall_interpolate_2d_vec
        write(*,*) "ncall_twister            ", ncall_twister      
        write(*,*) "ncall_twisterd           ", ncall_twisterd     
        write(*,*) "ncall_makeDmat           ", ncall_makeDmat    
        write(*,*) "ncall_rotate_fb          ", ncall_rotate_fb 

        write(*,*) "ncall_assemble_sVNL      ", ncall_assemble_sVNL 
        write(*,*) "ncall_assemble_2c        ", ncall_assemble_2c 
        write(*,*) "ncall_assemble_2c_PP     ", ncall_assemble_2c_PP 
        write(*,*) "ncall_average_ca_rho     ", ncall_average_ca_rho 
        write(*,*) "ncall_average_rho        ", ncall_average_rho 
        write(*,*) "ncall_assemble_snxc_on   ", ncall_assemble_snxc_on 
        write(*,*) "ncall_assemble_snxc_off  ", ncall_assemble_snxc_off 
        write(*,*) "ncall_assemble_olsxc_on  ", ncall_assemble_olsxc_on 
        write(*,*) "ncall_assemble_olsxc_off ", ncall_assemble_olsxc_off 
        write(*,*) "ncall_assemble_ca_2c     ", ncall_assemble_ca_2c 
        write(*,*) "ncall_assemble_3c        ", ncall_assemble_3c 
        write(*,*) "ncall_assemble_3c_PP     ", ncall_assemble_3c_PP 
        write(*,*) "ncall_assemble_ca_3c     ", ncall_assemble_ca_3c 
        write(*,*) "ncall_assemble_lr        ", ncall_assemble_lr 

        do i = 1,max_interaction
            write(*,*) "ncall_interpolate_1d_inter [",i,"] ", ncall_interpolate_1d_inter(i)
        enddo
        do i = 1,4
            write(*,*) "ncall_interpolate_2d_inter [",i,"] ", ncall_interpolate_2d_inter(i)
        enddo
    end subroutine write_ncall

    subroutine writeclean_ncall()
        call write_ncall()
        call clean_ncall()
    end subroutine writeclean_ncall

    subroutine timer_init()
        timer_index = 1
        timer_starts(:) = -1
        timer_sums  (:) = 0
    end subroutine timer_init

    subroutine timer_start_i(i)
        integer,                 intent (in) :: i
        if( timer_starts(i) .gt. 0 ) then
            write(*,*) "timer_start_i(",i,") ERORROR timer slot is occupied"
            stop
        end if 
        call cpu_time( timer_starts(i) )
    end subroutine timer_start_i

    subroutine timer_stop_i(i)
        integer,                 intent (in) :: i
        real :: tmp
        call cpu_time( tmp )
        timer_sums(i) = timer_sums(i) + tmp - timer_starts(i)
        timer_starts(i) = -1
    end subroutine timer_stop_i

    !subroutine timer_start()
    !    if( timer_index .gt. timer_index_max) then
    !        write(*,*) "timer_start() : timer_index > timer_index_max ", timer_index, timer_index_max
    !        stop
    !    end if 
    !    call cpu_time( timer_starts(timer_index) )
    !    timer_index = timer_index + 1
    !end subroutine timer_start

    !subroutine timer_stop()
    !    if( timer_index .lt. 11) then
    !        write(*,*) "timer_end() : timer_index < 1 ", timer_index
    !        stop
    !    end if 
    !    call cpu_time( timer_starts(timer_index) )
    !    timer_index = timer_index - 1
    !end subroutine timer_stop

end module timing