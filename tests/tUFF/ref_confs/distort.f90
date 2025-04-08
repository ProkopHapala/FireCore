program distort
  implicit none
  real(8), parameter :: tol = 1.d-9
  integer :: i, j, k, l, atom_type
  integer :: natoms, nbonds, nangles, ndihedrals, nimpropers
  integer :: natom_types, nbond_types, nangle_types, ndihedral_types, nimproper_types
  integer :: bonds(3), angles(4), dihedrals(5), impropers(5)
  real(8) :: xlo, xhi, ylo, yhi, zlo, zhi
  real(8) :: dist_max, m, q, pos(3), r, coeffs(5)
  character(200) :: file_in, file_data, file_xyz, comment, comment2
  character(3), allocatable :: atom_labels(:)
  
  read(*,*) file_in, file_data, file_xyz, dist_max
  
  open(unit=10,file=file_in)
  open(unit=100,file=file_data)
  open(unit=200,file=file_xyz)

  ! preamble
  read(10,*)
  write(100,*)
  read(10,*)
  write(100,*)

  ! number of atoms, bonds, etc.
  read(10,*) natoms
  write(100,*) natoms, ' atoms'
  write(200,*) natoms
  read(10,*) natom_types
  write(100,*) natom_types, ' atom types'
  read(10,*) nbonds
  write(100,*) nbonds, ' bonds'
  read(10,*) nbond_types
  write(100,*) nbond_types, ' bond types'
  read(10,*) nangles
  write(100,*) nangles, ' angles'
  read(10,*) nangle_types
  write(100,*) nangle_types, ' angle types'
  read(10,*) ndihedrals
  write(100,*) ndihedrals, ' dihedrals'
  read(10,*) ndihedral_types
  write(100,*) ndihedral_types, ' dihedral types'
  read(10,*) nimpropers
  write(100,*) nimpropers, ' impropers'
  read(10,*) nimproper_types
  write(100,*) nimproper_types, ' improper types'

  ! box
  read(10,*)
  write(100,*)
  read(10,*) xlo, xhi
  write(100,*) xlo, xhi, ' xlo xhi'
  read(10,*) ylo, yhi
  write(100,*) ylo, yhi, ' ylo yhi'
  read(10,*) zlo, zhi
  write(100,*) zlo, zhi, ' zlo zhi'
  write(200,*) 'lvs ', xhi-xlo,0.d0,0.d0, 0.d0,yhi-ylo,0.d0, 0.d0,0.d0,zhi-zlo
  
  ! masses
  read(10,*)
  write(100,*) 
  read(10,*) comment
  if ( comment /= 'Masses' ) stop 'Masses'
  write(100,*) 'Masses'
  read(10,*)
  write(100,*)
  allocate ( atom_labels(natom_types) )
  do i = 1, natom_types
     read(10,*) j, m
     write(100,*) j, m
     if ( abs( m - 12.0107d0 ) < tol ) atom_labels(i) = 'C'
     if ( abs( m - 15.9994d0 ) < tol ) atom_labels(i) = 'O'
     if ( abs( m - 1.00794d0 ) < tol ) atom_labels(i) = 'H'
  end do
  
  ! bond coeffs
  read(10,*)
  write(100,*) 
  read(10,*) comment, comment2
  if ( comment /= 'Bond' .or. comment2 /= 'Coeffs' ) stop 'Bond Coeffs'
  write(100,*) 'Bond Coeffs'
  read(10,*)
  write(100,*) 
  do i = 1, nbond_types
     read(10,*) k, (coeffs(j),j=1,2)
     write(100,*) k, (coeffs(j),j=1,2)
  end do

  ! angle coeffs
  read(10,*)
  write(100,*) 
  read(10,*) comment, comment2
  if ( comment /= 'Angle' .or. comment2 /= 'Coeffs' ) stop 'Angle Coeffs'
  write(100,*) 'Angle Coeffs'
  read(10,*)
  write(100,*) 
  do i = 1, nangle_types
     read(10,*) k, (coeffs(j),j=1,3)
     write(100,*) k, coeffs(1), (nint(coeffs(j)),j=2,3)
  end do

  ! dihedral coeffs
  read(10,*)
  write(100,*) 
  read(10,*) comment, comment2
  if ( comment /= 'Dihedral' .or. comment2 /= 'Coeffs' ) stop 'Dihedral Coeffs'
  write(100,*) 'Dihedral Coeffs'
  read(10,*)
  write(100,*) 
  do i = 1, ndihedral_types
     read(10,*) k, (coeffs(j),j=1,3)
     write(100,*) k, coeffs(1), (nint(coeffs(j)),j=2,3)
  end do

  ! improper coeffs
  read(10,*)
  write(100,*) 
  read(10,*) comment, comment2
  if ( comment /= 'Improper' .or. comment2 /= 'Coeffs' ) stop 'Improper Coeffs'
  write(100,*) 'Improper Coeffs'
  read(10,*)
  write(100,*) 
  do i = 1, nimproper_types
     read(10,*) k, (coeffs(j),j=1,5)
     write(100,*) k, (coeffs(j),j=1,4), nint(coeffs(5))
  end do

  ! atoms
  read(10,*)
  write(100,*) 
  read(10,*) comment
  if ( comment /= 'Atoms' ) stop 'Atoms'
  write(100,*) 'Atoms'
  read(10,*)
  write(100,*) 
  do i = 1, natoms
     read(10,*) k, l, atom_type, q, (pos(j),j=1,3)
     CALL RANDOM_NUMBER(r)
     do j = 1, 3
        pos(j) = pos(j) + ( r - 0.5d0 ) * 2.d0 * dist_max
     end do
     write(100,*) k, l, atom_type, q, (pos(j),j=1,3), 0, 0, 0
     write(200,*) atom_labels(atom_type), (pos(j),j=1,3), q
  end do

  ! bonds
  read(10,*)
  write(100,*) 
  read(10,*) comment
  if ( comment /= 'Bonds' ) stop 'Bonds'
  write(100,*) 'Bonds'
  read(10,*)
  write(100,*) 
  do i = 1, nbonds
     read(10,*) k, (bonds(j),j=1,3)
     write(100,*) k, (bonds(j),j=1,3)
  end do

  ! angles
  read(10,*)
  write(100,*) 
  read(10,*) comment
  if ( comment /= 'Angles' ) stop 'Angles'
  write(100,*) 'Angles'
  read(10,*)
  write(100,*) 
  do i = 1, nangles
     read(10,*) k, (angles(j),j=1,4)
     write(100,*) k, (angles(j),j=1,4)
  end do
  
  ! dihedrals
  read(10,*)
  write(100,*) 
  read(10,*) comment
  if ( comment /= 'Dihedrals' ) stop 'Dihedrals'
  write(100,*) 'Dihedrals'
  read(10,*)
  write(100,*) 
  do i = 1, ndihedrals
     read(10,*) k, (dihedrals(j),j=1,5)
     write(100,*) k, (dihedrals(j),j=1,5)
  end do
  
  ! impropers
  read(10,*)
  write(100,*) 
  read(10,*) comment
  if ( comment /= 'Impropers' ) stop 'Impropers'
  write(100,*) 'Impropers'
  read(10,*)
  write(100,*) 
  do i = 1, nimpropers
     read(10,*) k, (impropers(j),j=1,5)
     write(100,*) k, (impropers(j),j=1,5)
  end do

  close(10)
  close(100)
  close(200)

  stop
end program distort
