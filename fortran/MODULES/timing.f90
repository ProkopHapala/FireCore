
module timing
    implicit none

    integer, parameter :: timer_index_max = 100
    integer timer_index

    real timer_starts(timer_index_max)
    real timer_sums  (timer_index_max)

    contains

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