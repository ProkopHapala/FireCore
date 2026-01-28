program minmaxdiff
  implicit none
  real(8), parameter :: tol = 1.d-2
  integer :: i
  real(8) :: dmin, dmax, Emin, Emax, a, d, E
  character(200) :: filename

  read(*,*) filename, dmin, dmax

  ! read and find min/max
  Emin = 1.d99
  Emax = -1.d99
  open(unit=10,file=filename)
  do
     read(10,*,iostat=i) a, d, E
     if ( i < 0 ) exit
     if ( d < dmin-tol .or. d > dmax+tol ) cycle
     if ( E < Emin ) Emin = E
     if ( E > Emax ) Emax = E
  end do
  close(10)

  if ( Emin > 0.d0 ) then
     write(0,*) "ERROR: Emin must be < 0, got Emin= ", Emin
     stop
  end if

  if ( Emax < 0.d0 ) then
     write(0,*) "ERROR: Emax must be > 0, got Emax= ", Emax
     stop
  end if
  
  if ( abs(Emin) > abs(Emax) ) then
     Emax = -Emin
  else
     Emin = -Emax
  end if
  
  write(*,*) Emin, Emax

  stop
end program minmaxdiff
