program shift
  implicit none
  real(8), parameter :: tol = 1.d-9
  integer :: i, j, k, n
  integer :: natoms, nbonds, nangles, ndihedrals, nimpropers
  integer :: natom_types, nbond_types, nangle_types, ndihedral_types, nimproper_types
  integer :: bonds(3), angles(4), dihedrals(5), impropers(5)
  real(8) :: xlo, xhi, ylo, yhi, zlo, zhi, m, r, coeffs(5), shift_par
  character(7) :: shift_key
  character(200) :: file_in, file_data, file_xyzmol, file_xyzsub, comment, comment2
  integer, allocatable :: mol(:), atom_type(:)
  real(8), allocatable :: q(:), pos(:,:)
  character(3), allocatable :: atom_labels(:)
  
  read(*,*) file_in, file_data, file_xyzmol, file_xyzsub, shift_key, shift_par
  
  open(unit=10,file=file_in)
  open(unit=100,file=file_data)

  ! preamble
  read(10,*)
  write(100,*)
  read(10,*)
  write(100,*)

  ! number of atoms, bonds, etc.
  read(10,*) natoms
  write(100,*) natoms, ' atoms'
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
     if ( abs( m - 14.0067d0 ) < tol ) atom_labels(i) = 'N'
     if ( abs( m - 15.9994d0 ) < tol ) atom_labels(i) = 'O'
     if ( abs( m - 1.00794d0 ) < tol ) atom_labels(i) = 'H'
     if ( abs( m - 22.98976928d0 ) < tol ) atom_labels(i) = 'Na'
     if ( abs( m - 35.45d0 ) < tol ) atom_labels(i) = 'Cl'
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
  allocate ( mol(natoms) )
  allocate ( atom_type(natoms) )
  allocate ( q(natoms) )
  allocate ( pos(natoms,3) )
  n = 0
  do i = 1, natoms
     read(10,*) k, mol(i), atom_type(i), q(i), (pos(i,j),j=1,3)
     if ( mol(i) == 1 ) then
        n = n + 1
        if ( shift_key == 'distort' ) then
           CALL RANDOM_NUMBER(r)
           do j = 1, 3
              pos(i,j) = pos(i,j) + ( r - 0.5d0 ) * 2.d0 * shift_par
           end do
        else if ( shift_key == 'xscan' ) then
           pos(i,1) = pos(i,1) + shift_par
        else if ( shift_key == 'yscan' ) then
           pos(i,2) = pos(i,2) + shift_par
        else if ( shift_key == 'zscan' ) then
           pos(i,3) = pos(i,3) + shift_par
        end if
     end if
     write(100,*) k, mol(i), atom_type(i), q(i), (pos(i,j),j=1,3), 0, 0, 0
  end do

  open(unit=200,file=file_xyzmol)
  open(unit=300,file=file_xyzsub)
  write(200,*) n
  write(200,*) 'lvs ', xhi-xlo,0.d0,0.d0, 0.d0,yhi-ylo,0.d0, 0.d0,0.d0,zhi-zlo
  write(300,*) natoms-n
  write(300,'(a,9f13.3)') 'lvs ', xhi-xlo,0.d0,0.d0, 0.d0,yhi-ylo,0.d0, 0.d0,0.d0,zhi-zlo
  do i = 1, natoms
     if ( mol(i) == 1 ) then
        write(200,*) atom_labels(atom_type(i)), (pos(i,j),j=1,3), q(i)
     else
        write(300,*) atom_labels(atom_type(i)), (pos(i,j),j=1,3), q(i)
     end if
  end do
  close(200)
  close(300)
  
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
end program shift
