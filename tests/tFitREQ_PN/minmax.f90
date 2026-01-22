program minmax
  implicit none
  real(8), parameter :: tol = 1.d-2
  integer :: i, j, na, nd, n
  real(8) :: dmin, dmax, Emin, Emax, aa, aa_old
  character(200) :: filein, fileout
  real(8), allocatable :: a(:,:), d(:,:), E(:,:)

  read(*,*) filein, fileout, dmin, dmax

  ! read and store
  open(unit=10,file=filein)
  read(10,*) aa_old
  na = 1
  n = 1
  do
     read(10,*,iostat=i) aa
     if ( i < 0 ) exit
     if ( abs(aa-aa_old) > tol ) then
        if ( na == 1 ) then
           nd = n
        else
           if ( n /= nd ) then
              write(0,*) 'STOP'
              write(0,*) 'na = ', na
              write(0,*) 'nd = ', nd
              write(0,*) 'n = ', n
              stop
           end if
        end if
        na = na + 1
        n = 0
     end if
     n = n + 1
     aa_old = aa
  end do
  rewind(10)
  allocate( a(na,nd) )
  allocate( d(na,nd) )
  allocate( E(na,nd) )
  do i = 1, na
     do j = 1, nd
        read(10,*) a(i,j), d(i,j), E(i,j)
     end do
  end do
  close(10)

  ! find min/max
  Emin = 1.d99
  Emax = -1.d99
  do i = 1, na
     do j = 1, nd
        if ( d(i,j) < dmin-tol .or. d(i,j) > dmax+tol ) cycle
        if ( E(i,j) < Emin ) Emin = E(i,j)
        if ( E(i,j) > Emax ) Emax = E(i,j)
     end do
  end do
  write(*,*) Emin, Emax

  ! write for cart
  open(unit=11,file=fileout)
  do j = nd, 1, -1
     do i = 1, na
        if ( d(i,j) < dmin-tol .or. d(i,j) > dmax+tol ) cycle
        write(11,*) a(i,j), d(i,j), E(i,j)
     end do
     write(11,*)
  end do
  close(11)

  stop
end program minmax
