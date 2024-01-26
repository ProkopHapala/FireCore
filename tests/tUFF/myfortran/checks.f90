subroutine checks
  use conf, only: tol, natoms, ufftypes, neighbors, connectivity, bond_orders, &
       nbonds, bond_orders_int, set_atom, set_bond
  implicit none
  integer :: ia, j, ja, ib
  logical :: found

  ! sanity checks
  do ia = 1, natoms
     if ( .not. set_atom(ia) ) then
        write(0,*) 'ERROR CHECKS: atom ', ia, ' is not set'
        write(0,*) 'STOP'
        stop
     end if
     if ( ufftypes(ia) == '' ) then
        write(0,*) 'ERROR CHECKS: type for atom ', ia, ' is not set'
        write(0,*) 'STOP'
        stop
     end if
  end do
  do ib = 1, nbonds
     if ( .not. set_bond(ib) ) then
        write(0,*) 'ERROR CHECKS: bond ', ib, ' is not set'
        write(0,*) 'STOP'
        stop
     end if
     if ( bond_orders(ib) < 0.d0 ) then
        write(0,*) 'ERROR CHECKS: order for bond ', ib, ' is not set'
        write(0,*) 'STOP'
        stop
     end if
     if ( bond_orders_int(ib) < 0 ) then
        write(0,*) 'ERROR CHECKS: integer order for bond ', ib, ' is not set'
        write(0,*) 'STOP'
        stop
     end if
  end do
  
  do ia = 1, natoms
     if ( ufftypes(ia)(3:3) == 'R' ) then
        ! check that resonant atoms have at least one 1.5 bond
        found = .false.
        do j = 1, neighbors(ia)
           ja = connectivity(ia,j)
           call find_bond ( ia, ja, ib )
           if ( abs(bond_orders(ib)-1.5d0) < tol ) then
              found = .true.
              exit
           end if
        end do
        if ( .not. found ) then
           write(0,*) 'ERROR CHECKS: atom ', ia, ' is resonant but it has no 1.5 bonds'
           do ja = 1, neighbors(ia)
              call find_bond ( ia, connectivity(ia,ja), ib )
              write(0,*) '  bond ', ib, ' atoms ', ufftypes(ia), ia, ufftypes(connectivity(ia,ja)), connectivity(ia,ja), &
                   ' bond order ', bond_orders(ib), bond_orders_int(ib)
           end do
           write(0,*) 'STOP'
           stop
        end if
     else if ( ufftypes(ia)(3:3) == '2' ) then
        ! check that sp2 atoms have one localized double bond
        found = .false.
        do j = 1, neighbors(ia)
           ja = connectivity(ia,j)
           call find_bond ( ia, ja, ib )
           if ( abs(bond_orders(ib)-2.d0) < tol ) then
              if ( found ) then
                 write(0,*) 'ERROR CHECKS: atom ', ia, ' is sp2 but it has more than one double bond'
                 do ja = 1, neighbors(ia)
                    call find_bond ( ia, connectivity(ia,ja), ib )
                    write(0,*) '  bond ', ib, ' atoms ', ufftypes(ia), ia, ufftypes(connectivity(ia,ja)), connectivity(ia,ja), &
                         ' bond order ', bond_orders(ib), bond_orders_int(ib)
                 end do
                 write(0,*) 'STOP'
                 stop
              end if
              found = .true.
           end if
        end do
        if ( .not. found ) then
           found = .false.
           do j = 1, neighbors(ia)
              ja = connectivity(ia,j)
              call find_bond ( ia, ja, ib )
              if ( abs(bond_orders(ib)-1.5d0) < tol ) found = .true.
           end do
           if ( found ) then
              write(*,*) 'WARNING CHECKS: atom ', ia, ' is sp2 and it has no double bonds, only 1.5 bonds'
           else
              write(0,*) 'ERROR CHECKS: atom ', ia, ' is sp2 but it has no double bonds'
              do ja = 1, neighbors(ia)
                 call find_bond ( ia, connectivity(ia,ja), ib )
                 write(0,*) '  bond ', ib, ' atoms ', ufftypes(ia), ia, ufftypes(connectivity(ia,ja)), connectivity(ia,ja), &
                      ' bond order ', bond_orders(ib), bond_orders_int(ib)
              end do
              write(0,*) 'STOP'
              stop
           end if
        end if
     end if
  end do
        
  ! check that 1.5 bonds must be either 1 or 2 in the limit resonance structure
  do ib = 1, nbonds
     if ( abs(bond_orders(ib)-1.5d0) < tol ) then
        if ( bond_orders_int(ib) /= 1 .and. bond_orders_int(ib) /= 2 ) then
           write(0,*) 'ERROR CHECKS: bond ', ib, ' is 1.5 but in the limit resonance structure is neither single nor double'
           do ja = 1, neighbors(ia)
              call find_bond ( ia, connectivity(ia,ja), ib )
              write(0,*) '  bond ', ib, ' atoms ', ufftypes(ia), ia, ufftypes(connectivity(ia,ja)), connectivity(ia,ja), &
                   ' bond order ', bond_orders(ib), bond_orders_int(ib)
           end do
           write(0,*) 'STOP'
           stop
        end if
     end if
  end do
  
  return
end subroutine checks
