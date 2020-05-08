


module workmat

    integer lwork, lrwork, liwork
    complex, allocatable, dimension (:) :: work
    real,    allocatable, dimension (:) :: rwork
    integer, allocatable, dimension (:) :: iwork

    real, allocatable, dimension (:) :: slam

    complex, allocatable, dimension (:, :)    :: xxxx
    complex, allocatable, dimension (:, :, :) :: sm12_save

end module workmat
