subroutine fix_saturation
  use conf, only: tol, natoms, neighbors, connectivity, valence_target, bond_orders, set_bond, tipo, bond_orders_int
  implicit none
  integer :: ia, ib, na, j, ja, nset
  real(8) :: val
  logical :: changed, ok
  
  do
     changed = .false.
     do ia = 1, natoms
        na = neighbors(ia)
        ! if the atom has all bonds set, skip
        ok = .true.
        do j = 1, na
           ja = connectivity(ia,j)
           call find_bond ( ia, ja, ib )
           if ( bond_orders(ib) < 0.d0 ) ok = .false.
        end do
        if ( ok ) cycle
        ! compute number of already set bonds and corresponding atom valence
        nset = 0
        val = 0.d0
        do j = 1, na
           ja = connectivity(ia,j)
           call find_bond ( ia, ja, ib )
           if ( bond_orders(ib) > 0.d0 ) then
              nset = nset + 1
              !val = val + bond_orders(ib)
              val = val + dble(bond_orders_int(ib))
           end if
        end do
        val = dble(valence_target(ia)) - val
        ! only one bond to be set
        if ( nset == na-1 ) then
           ! check for weird valence
           if ( abs( val - dble(nint(val)) ) > tol .or. abs(val) < tol ) then
              write(0,*) 'ERROR CHECKSATURATION: atom ', tipo(ia), ia, ' would have a valence of ', val
              write(0,*) 'STOP'
              stop
           end if
           changed = .true.
           do j = 1, na
              ja = connectivity(ia,j)
              call find_bond ( ia, ja, ib )
              if ( bond_orders(ib) < 0.d0 ) then
                 bond_orders(ib) = val
                 set_bond(ib) = .true.
                 exit
              end if
           end do
        ! or multiple bonds to be set, but they are all single
        elseif ( nint(val) == na-nset ) then
           changed = .true.
           do j = 1, na
              ja = connectivity(ia,j)
              call find_bond ( ia, ja, ib )
              if ( bond_orders(ib) < 0.d0 ) then
                 bond_orders(ib) = 1.d0
                 set_bond(ib) = .true.
              end if
           end do
        end if
     end do
     if ( .not. changed ) return
  end do
 
  return
end subroutine fix_saturation
