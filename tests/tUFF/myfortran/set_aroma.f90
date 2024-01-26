subroutine set_aroma ( natoms, atoms )
  use conf, only: tol, bond_orders, neighbors, connectivity, tipo, &
       set_bond, ufftypes, set_atom
  implicit none
  integer, intent(in) :: natoms, atoms(natoms+1)
  integer :: i, ia, ib, j, jj
  logical :: ok
  
  ! set atoms
  do ia = 1, natoms
     if ( set_atom(atoms(ia)) ) then
        if ( ufftypes(atoms(ia))(3:3) /= 'R' ) write(*,*) 'WARNING SETAROMA: atom ', tipo(atoms(ia)), atoms(ia), &
             ' would be set to resonant but it has already a type of ', ufftypes(atoms(ia))
        cycle
     end if
     if ( tipo(atoms(ia)) == 'C' ) then
        ufftypes(atoms(ia)) = 'C_R'
        set_atom(atoms(ia)) = .true.
     end if
     if ( tipo(atoms(ia)) == 'N' ) then
        ufftypes(atoms(ia)) = 'N_R'
        set_atom(atoms(ia)) = .true.
     end if
     if ( tipo(atoms(ia)) == 'O' ) then
        ufftypes(atoms(ia)) = 'O_R'
        set_atom(atoms(ia)) = .true.
     end if
  end do

  ! set ring bonds
  do ia = 1, natoms
     call find_bond ( atoms(ia), atoms(ia+1), ib )
     if ( set_bond(ib) ) then
        if ( abs(bond_orders(ib)-1.5d0) > tol ) write(*,*) 'WARNING SETAROMA: bond ', ib, ' between atoms ', &
             (tipo(atoms(ia+j)),atoms(ia+j),j=0,1), ' would be set to 1.5 but it has already a bond order of ', bond_orders(ib)
        cycle
     end if
     bond_orders(ib) = 1.5d0
     set_bond(ib) = .true.
  end do

  ! set other bonds
  do ia = 1, natoms
     do j = 1, neighbors(atoms(ia))
        jj = connectivity(atoms(ia),j)
        ok = .true.
        do i = 1, natoms
           if ( jj == atoms(i) ) ok = .false.
        end do
        if ( .not. ok ) cycle
        call find_bond ( atoms(ia), jj, ib )
        if ( tipo(jj) == 'O' .and. neighbors(jj) == 1 ) then ! delocalized carbonyl
           if ( set_atom(jj) ) then
              if ( ufftypes(jj)(3:3) /= 'R' ) write(*,*) 'WARNING SETAROMA: carbonyl atom ', tipo(jj), jj, &
                   ' would be set to resonant but it has already a type of ', ufftypes(jj)
              cycle
           end if
           if ( set_bond(ib) ) then
              if ( abs(bond_orders(ib)-1.5d0) > tol ) write(*,*) 'WARNING SETAROMA: lateral bond ', ib, ' between atoms ', &
                   tipo(atoms(ia)),atoms(ia),tipo(jj),jj, ' would be set to 1.5 but it has already a bond order of ', bond_orders(ib)
              cycle
           end if
           ufftypes(jj) = 'O_R'
           set_atom(jj) = .true.
           bond_orders(ib) = 1.5d0
           set_bond(ib) = .true.
!        else
!           bond_orders(ib) = 1.d0
!           set_bond(ib) = .true.
        end if
     end do
  end do
    
  return
end subroutine set_aroma
