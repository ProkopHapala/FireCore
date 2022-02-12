
module timing
    implicit none

    integer, parameter :: timer_index_max = 100
    integer timer_index

    integer ncall_doscentros
    integer ncall_doscentrosS
    integer ncall_trescentros
    integer ncall_trescentrosS
    integer ncall_epsilon
    integer ncall_deps2cent
    integer ncall_interpolate_1d
    integer ncall_interpolate_2d
    integer ncall_twister
    integer ncall_twisterd
    integer ncall_makeDmat
    integer ncall_rotate_fb

    real timer_starts(timer_index_max)
    real timer_sums  (timer_index_max)

    contains

    subroutine clean_ncall()
        ncall_doscentros   =0
        ncall_doscentrosS  =0
        ncall_trescentros  =0
        ncall_trescentrosS =0
        ncall_epsilon      =0
        ncall_deps2cent    =0
        ncall_interpolate_1d =0
        ncall_interpolate_2d =0
        ncall_twister      =0
        ncall_twisterd     =0
        ncall_makeDmat     =0
        ncall_rotate_fb    =0
    end subroutine clean_ncall

    subroutine write_ncall()
        write(*,*) "ncall_doscentros",     ncall_doscentros
        write(*,*) "ncall_doscentrosS",    ncall_doscentrosS  
        write(*,*) "ncall_trescentros",    ncall_trescentros  
        write(*,*) "ncall_trescentrosS",   ncall_trescentrosS 
        write(*,*) "ncall_epsilon",        ncall_epsilon      
        write(*,*) "ncall_deps2cent",      ncall_deps2cent    
        write(*,*) "ncall_interpolate_1d", ncall_interpolate_1d 
        write(*,*) "ncall_interpolate_2d", ncall_interpolate_2d 
        write(*,*) "ncall_twister",        ncall_twister      
        write(*,*) "ncall_twisterd",       ncall_twisterd     
        write(*,*) "ncall_makeDmat",       ncall_makeDmat    
        write(*,*) "ncall_rotate_fb",      ncall_rotate_fb 
    end subroutine write_ncall

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