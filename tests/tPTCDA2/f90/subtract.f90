program subtract
  implicit none
  real(8), parameter :: tol = 1.d-8
  integer :: i, j, k, natoms
  real(8) :: e_current, e_offset
  real(8), allocatable :: pos_current(:,:), f_current(:,:), pos_offset(:,:), f_offset(:,:)
  
  read(*,*) natoms

  allocate(pos_current(natoms,3))
  allocate(f_current(natoms,3))
  open(unit=10,file='f_current.txt')
  do i = 1, natoms
     read(10,*) k, (pos_current(i,j),j=1,3), (f_current(i,j),j=1,3)
     if ( k /= i ) stop 'superscio'
  end do
  close(10)
  open(unit=10,file='e_current.txt')
  read(10,*) e_current
  close(10)

  allocate(pos_offset(natoms,3))
  allocate(f_offset(natoms,3))
  open(unit=10,file='f_offset.txt')
  do i = 1, natoms
     read(10,*) k, (pos_offset(i,j),j=1,3), (f_offset(i,j),j=1,3)
     if ( k /= i ) stop 'superscio'
  end do
  close(10)
  open(unit=10,file='e_offset.txt')
  read(10,*) e_offset
  close(10)

  do i = 1, natoms
     do j = 1, 3
        if ( abs( pos_current(i,j) - pos_offset(i,j) ) > tol ) stop 'jamme'
     end do
  end do

  open(unit=10,file='f_lammps.txt')
  do i = 1, natoms
     write(10,*) i, (pos_current(i,j),j=1,3), (f_current(i,j)-f_offset(i,j),j=1,3)
  end do
  close(10)
  
  open(unit=10,file='e_lammps.txt')
  write(10,*) e_current-e_offset
  close(10)

  stop
end program subtract
