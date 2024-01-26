subroutine write_data ( output, output_qeq, tagcharge )
  use constants, only: mass_H, mass_C, mass_N, mass_O
  use conf, only: ufftypes, natoms, nbonds, bonds, nangles, angles, ntorsions, torsions, ninversions, inversions, side, pos, q, &
       natomtypes, nbondtypes, nangletypes, ntorsiontypes, ninversiontypes, atomtypesmap, bondtypesmap, angletypesmap, torsiontypesmap, inversiontypesmap
  use pars, only: vdw_pars, bond_pars, angle_tags, angle_pars, torsion_pars, inversion_pars
  use uff, only: chi, J_I, eta, Zstar_I, chi_reax, eta_reax, gamma_reax
  implicit none
  character(20), intent(in) :: tagcharge
  character(200), intent(in) :: output, output_qeq
  integer :: i, j, k
  
  open(unit=100,file=output)
  write(100,*)
  write(100,*)
  write(100,*) natoms, ' atoms'
  write(100,*) natomtypes, ' atom types'
  write(100,*) nbonds, ' bonds'
  write(100,*) nbondtypes, ' bond types'
  write(100,*) nangles, ' angles'
  write(100,*) nangletypes, ' angle types'
  write(100,*) ntorsions, ' dihedrals'
  write(100,*) ntorsiontypes, ' dihedral types'
  write(100,*) ninversions, ' impropers'
  write(100,*) ninversiontypes, ' improper types'
  write(100,*)
  write(100,*) 0.d0, side(1,1), ' xlo xhi'
  write(100,*) 0.d0, side(2,2), ' ylo yhi'
  write(100,*) 0.d0, side(3,3), ' zlo zhi'
  write(100,*) side(2,1), side(3,1), side(3,2), ' xy xz yz'
  write(100,*)
  write(100,*) 'Masses'
  write(100,*)
  do i = 1, natomtypes
     do k = 1, natoms
        if ( atomtypesmap(k) == i ) exit
     end do
     if ( k > natoms ) stop 'superscio'
     if ( ufftypes(k)(1:1) == 'H' ) then
        write(100,*) i, mass_H, ' # ', ufftypes(k)
     else if ( ufftypes(k)(1:1) == 'C' ) then
        write(100,*) i, mass_C, ' # ', ufftypes(k)
     else if ( ufftypes(k)(1:1) == 'N' ) then
        write(100,*) i, mass_N, ' # ', ufftypes(k)
     else if ( ufftypes(k)(1:1) == 'O' ) then
        write(100,*) i, mass_O, ' # ', ufftypes(k)
     end if
  end do
  if ( nbondtypes > 0 ) then
     write(100,*)
     write(100,*) 'Bond Coeffs'
     write(100,*)
     do i = 1, nbondtypes
        do k = 1, nbonds
           if ( bondtypesmap(k) == i ) exit
        end do
        if ( k > nbonds ) stop 'superscio'
        write(100,*) i, (bond_pars(i,j),j=1,2), ' # ', (ufftypes(bonds(k,j)),j=1,2)
     end do
  end if
  if ( nangletypes > 0 ) then
     write(100,*)
     write(100,*) 'Angle Coeffs'
     write(100,*)
     do i = 1, nangletypes
        do k = 1, nangles
           if ( angletypesmap(k) == i ) exit
        end do
        if ( k > nangles ) stop 'superscio'
        if ( angle_tags(i) == 'cosine/periodic' ) then
           write(100,*) i, angle_tags(i), angle_pars(i,1), nint(angle_pars(i,2)), nint(angle_pars(i,3)), ' # ', (ufftypes(angles(k,j)),j=1,3)
        else if ( angle_tags(i) == 'fourier' ) then
           write(100,*) i, angle_tags(i), (angle_pars(i,j),j=1,4), ' # ', (ufftypes(angles(k,j)),j=1,3)
        end if
     end do
  end if
  if ( ntorsiontypes > 0 ) then
     write(100,*)
     write(100,*) 'Dihedral Coeffs'
     write(100,*)
     do i = 1, ntorsiontypes
        do k = 1, ntorsions
           if ( torsiontypesmap(k) == i ) exit
        end do
        if ( k > ntorsions ) stop 'superscio'
        write(100,*) i, torsion_pars(i,1), nint(torsion_pars(i,2)), nint(torsion_pars(i,3)), ' # ', (ufftypes(torsions(k,j)),j=1,4)
     end do
  end if
  if ( ninversiontypes > 0 ) then
     write(100,*)
     write(100,*) 'Improper Coeffs'
     write(100,*)
     do i = 1, ninversiontypes
        do k = 1, ninversions
           if ( inversiontypesmap(k) == i ) exit
        end do
        if ( k > ninversions ) stop 'superscio'
        write(100,*) i, (inversion_pars(i,j),j=1,4), 0, ' # ', (ufftypes(inversions(k,j)),j=1,4)
     end do
  end if
  write(100,*)
  write(100,*) 'Pair Coeffs'
  write(100,*)
  do i = 1, natomtypes
     do k = 1, natoms
        if ( atomtypesmap(k) == i ) exit
     end do
     if ( k > natoms ) stop 'superscio'
     write(100,*) i, (vdw_pars(i,j),j=1,2), ' # ', ufftypes(k), ufftypes(k)
  end do
  write(100,*)
  write(100,*) 'Atoms # full'
  write(100,*)
  do i = 1, natoms
     if ( tagcharge == 'nocharges' ) then
        write(100,*) i, 1, atomtypesmap(i), 0.d0, (pos(i,j),j=1,3), 0, 0, 0
     else
        write(100,*) i, 1, atomtypesmap(i), q(i), (pos(i,j),j=1,3), 0, 0, 0
     end if
  end do
  if ( nbonds > 0 ) then
     write(100,*)
     write(100,*) 'Bonds'
     write(100,*)
     do i = 1, nbonds
        write(100,*) i, bondtypesmap(i), (bonds(i,j),j=1,2)
     end do
  end if
  if ( nangles > 0 ) then
     write(100,*)
     write(100,*) 'Angles'
     write(100,*)
     do i = 1, nangles
        write(100,*) i, angletypesmap(i), (angles(i,j),j=1,3)
     end do
  end if
  if ( ntorsions > 0 ) then
     write(100,*)
     write(100,*) 'Dihedrals'
     write(100,*)
     do i = 1, ntorsions
        write(100,*) i, torsiontypesmap(i), (torsions(i,j),j=1,4)
     end do
  end if
  if ( ninversions > 0 ) then
     write(100,*)
     write(100,*) 'Impropers'
     write(100,*)
     do i = 1, ninversions
        write(100,*) i, inversiontypesmap(i), (inversions(i,j),j=1,4)
     end do
  end if
  close(100)

  if ( tagcharge /= 'nocharges' .and. tagcharge /= 'xyzcharges' ) then
     open(unit=110,file=output_qeq)
     !write(110,'(a)') '# type chi(kcal/mol) eta(kcal/mol) gamma(1/Ang) zeta(1/Ang) qcore(e)'
     do i = 1, natomtypes
        do k = 1, natoms
           if ( atomtypesmap(k) == i ) exit
        end do
        if ( k > natoms ) stop 'superscio'
        call find_I ( ufftypes(k), j )
        if ( tagcharge == 'qeq-point' .or. tagcharge == 'qeq-slater' ) then
           write(110,*) i, chi(j), J_I(j), 0.d0, eta(j), Zstar_I(j)
        else if ( tagcharge == 'qeq-shielded' ) then
           write(110,*) i, chi(j), J_I(j), gamma_reax(j), 0.d0, 0.d0
        else if ( tagcharge == 'qeq-reaxff' ) then
           write(110,*) i, chi_reax(j), eta_reax(j)*2.d0, gamma_reax(j)
        else
           write(0,*) 'wrong charge keyword'
           stop
        end if
     end do
     close(110)
  end if
  
  open(unit=200,file='elements.tmp')
  do i = 1, natomtypes
     do k = 1, natoms
        if ( atomtypesmap(k) == i ) exit
     end do
     if ( k > natoms ) stop 'superscio'
     write(200,'(a2)',advance='no') ufftypes(k)(1:1)//' '
  end do
  write(200,*) ''
  close(200)
  
  return
end subroutine write_data
