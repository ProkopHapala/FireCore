program calc_rmse
  implicit none
  integer :: i, n
  real(8) :: a, d, E, rmse
  character(200) :: filename

  read(*,*) filename

  ! read and find min/max
  rmse = 0.d0
  n = 0
  open(unit=10,file=filename)
  do
     read(10,*,iostat=i) a, d, E
     if ( i < 0 ) exit
     rmse = rmse + E * E
     n = n + 1
  end do
  close(10)
  
  write(*,'(f6.2)') sqrt( rmse / dble(n) )

  stop
end program calc_rmse
