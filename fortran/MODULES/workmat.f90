


module workmat

    integer lwork, lrwork, liwork
    complex, allocatable, dimension (:) :: work
    real,    allocatable, dimension (:) :: rwork
    integer, allocatable, dimension (:) :: iwork

    real, allocatable, dimension (:) :: slam

    complex, allocatable, dimension (:, :)    :: xxxx, zzzz
    complex, allocatable, dimension (:, :, :) :: sm12_save

contains 

    subroutine alloc_work( norbitals, divide )
        ! ===== Parameters
        integer,intent (in)  :: norbitals
        logical, intent(in) :: divide  ! do we divide and conquer?
        ! ===== Variables
        ! ===== Body
        !write(*,*) "DEBUG alloc_work norbitals, divide ", norbitals, divide 
        allocate ( slam(norbitals) )
        allocate ( xxxx(norbitals,norbitals) )
        allocate ( zzzz(norbitals,norbitals) )
        liwork = 10*norbitals
        allocate (iwork(liwork))
        if (divide) then
            lwork  = 100*norbitals +   norbitals*norbitals
            lrwork = 100*norbitals + 3*norbitals*norbitals ! for old versions of cheevd
            !liwork = 10*norbitals
            allocate (work(lwork))
            allocate (rwork(lrwork))
        else
        !! original
        !   lwork = norbitals*norbitals ! Use xxxx, yyyy and zzzz for work area
        !   lrwork = 3*norbitals
        !   allocate (rwork(lrwork))
            lwork = 1
            allocate (work(lwork))
            lrwork = 3*norbitals - 2
            allocate (rwork(lrwork))
        end if
    end subroutine alloc_work
        
    subroutine dealloc_work( )
        ! ===== Parameters
        !integer ,intent (in) :: divide
        ! ===== Variables
        ! ===== Body
        deallocate ( work )
        deallocate (iwork )
        deallocate (rwork )
        deallocate ( xxxx )
        deallocate ( zzzz )
        deallocate ( slam )
    end subroutine dealloc_work
        
end module workmat
