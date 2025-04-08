program compare
  implicit none
  real(8), parameter :: tol = 1.d-8
  integer :: i, j, k, natoms
  real(8) :: e_lammps, e_firecore
  real(8), allocatable :: pos_lammps(:,:), f_lammps(:,:), pos_firecore(:,:), f_firecore(:,:)
  
  read(*,*) natoms

  allocate(pos_lammps(natoms,3))
  allocate(f_lammps(natoms,3))
  open(unit=10,file='f_lammps.txt')
  do i = 1, natoms
     read(10,*) k, (pos_lammps(i,j),j=1,3), (f_lammps(i,j),j=1,3)
     if ( k /= i ) stop 'superscio'
  end do
  close(10)
  open(unit=10,file='e_lammps.txt')
  read(10,*) e_lammps
  close(10)

  allocate(pos_firecore(natoms,3))
  allocate(f_firecore(natoms,3))
  open(unit=10,file='f_firecore.txt')
  do i = 1, natoms
     read(10,*) k, (pos_firecore(i,j),j=1,3), (f_firecore(i,j),j=1,3)
     if ( k /= i ) stop 'superscio'
  end do
  close(10)
  open(unit=10,file='e_firecore.txt')
  read(10,*) e_firecore
  close(10)

  do i = 1, natoms
     do j = 1, 3
        if ( abs( pos_lammps(i,j) - pos_firecore(i,j) ) > tol ) then
           write(*,*) 'ATOM POSITIONS DO NOT MATCH'
           write(*,*) 'ATOM ', i, ' COMPONENT ', j
           write(*,*) 'POS LAMMPS   = ', pos_lammps(i,j)
           write(*,*) 'POS FIRECORE = ', pos_firecore(i,j)
           stop
        end if
     end do
  end do

  do i = 1, natoms
     do j = 1, 3
        if ( abs(f_lammps(i,j)) > tiny(1.d0) ) then
           write(*,*) 'FORCE ', i, j, f_lammps(i,j), f_firecore(i,j), abs((f_lammps(i,j)-f_firecore(i,j))/f_lammps(i,j))
        else
           write(*,*) 'FORCE ', i, j, f_lammps(i,j), f_firecore(i,j), abs(f_firecore(i,j))
        end if
     end do
  end do

  if ( abs(e_lammps) > tiny(1.d0) ) then
        write(*,*) 'ENERGY ', e_lammps, e_firecore, abs((e_lammps-e_firecore)/e_lammps)
  else
        write(*,*) 'ENERGY ', e_lammps, e_firecore, abs(e_firecore)
  end if

  stop
end program compare
