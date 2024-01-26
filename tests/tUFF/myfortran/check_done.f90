subroutine check_done ( bond_orders, done )
  use conf, only: nbonds
  implicit none
  integer, intent(in) :: bond_orders(nbonds)
  logical, intent(out) :: done
  integer :: ib

  done = .true.

  ! check if the bond orders are all set
  do ib = 1, nbonds
     if ( bond_orders(ib) < 0 ) then
        done = .false.
        return
     end if
  end do
  
  return
end subroutine check_done

