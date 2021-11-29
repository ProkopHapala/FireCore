
module debug
    implicit none
    contains

subroutine debug_writeArray_1i( name, arr, n )
    ! ===== Parameters
    character(len=*)      , intent(in)  :: name
    integer, dimension( : ), intent (in) :: arr
    integer,                 intent (in) :: n
    ! ===== Variables
    integer i
    ! ===== Body
    do i=1,n
        write(*,'(A,i6,A,i6)') "name[",i,"]: ", arr(i)
    end do
end subroutine debug_writeArray_1i

end module debug