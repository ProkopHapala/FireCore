subroutine check_all ( bond_orders, done )
  use conf, only: nbonds
  implicit none
  integer, intent(inout) :: bond_orders(nbonds)
  logical, intent(out) :: done
  logical :: changed

  ! check saturation
  do 
     changed = .false.
     call check_sat ( bond_orders, changed )
     if ( .not. changed ) exit
  end do
  
  ! check if we are done
  call check_done ( bond_orders, done )

  return
end subroutine check_all
