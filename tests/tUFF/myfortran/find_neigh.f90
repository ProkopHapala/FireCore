subroutine find_neigh
  use conf, only: natoms, nbonds, bonds, tipo, neighbors, connectivity, side, pos, &
       valence_target, bond_orders, bond_orders_int, set_bond!, adhoc_bond
  implicit none
  integer :: i, j, bonds_tmp(natoms*4,2)
  real(8) :: d, deq, side_inv(3,3), pos_frac(natoms,3)

  ! get fractional coordinates
  call inverse ( side, side_inv )
  do i = 1, natoms
     do j = 1, 3
        pos_frac(i,j) = side_inv(1,j)*pos(i,1) + side_inv(2,j)*pos(i,2) + side_inv(3,j)*pos(i,3)
     end do
  end do

  ! find bonds
  nbonds = 0
  do i = 1, natoms-1
     do j = i+1, natoms
        if ( tipo(i) == 'H' .and. tipo(j) == 'H' ) cycle ! exception to avoid H-H bonds
        call calc_dist ( pos_frac(i,:), pos_frac(j,:), d )
        call calc_dist_eq ( i, j, deq )
        if ( d < deq ) then
           nbonds = nbonds + 1
           bonds_tmp(nbonds,1) = i
           bonds_tmp(nbonds,2) = j
           neighbors(i) = neighbors(i) + 1
           if ( neighbors(i) > valence_target(i) ) then
              write(0,*) 'ERROR FINDNEIGH weird topology: atom ', tipo(i), i, ' has too many  neighbors!'
              write(0,*) 'ERROR FINDNEIGH valence: ', valence_target(i)
              write(0,*) 'ERROR FINDNEIGH neighbors: ', connectivity(i,:), j
              write(0,*) 'STOP'
              stop
           end if
           connectivity(i,neighbors(i)) = j
           neighbors(j) = neighbors(j) + 1
           if ( neighbors(j) > valence_target(j) ) then
              write(0,*) 'ERROR FINDNEIGH weird topology: atom ', tipo(j), j, ' has too many  neighbors!'
              write(0,*) 'ERROR FINDNEIGH valence: ', valence_target(j)
              write(0,*) 'ERROR FINDNEIGH neighbors: ', connectivity(j,:), i
              write(0,*) 'STOP'
              stop
           end if
           connectivity(j,neighbors(j)) = i
        end if
     end do
  end do

  ! allocate and initialize/store
  allocate ( bonds(nbonds,2), stat = i )
  if ( i /= 0 ) stop 'allocation BONDS'
  do i = 1, nbonds
     do j = 1, 2
        bonds(i,j) = bonds_tmp(i,j)
     end do
  end do
  allocate ( bond_orders(nbonds), stat = i )
  if ( i /= 0 ) stop 'allocation BOND_ORDERS'
  bond_orders(:) = -1.d0
  allocate ( bond_orders_int(nbonds), stat = i )
  if ( i /= 0 ) stop 'allocation BOND_ORDERS_INT'
  bond_orders_int(:) = -1
  allocate ( set_bond(nbonds), stat = i )
  if ( i /= 0 ) stop 'allocation SET_BOND'
  set_bond(:) = .false.
!  allocate ( adhoc_bond(nbonds), stat = i )
!  if ( i /= 0 ) stop 'allocation ADHOC_BOND'
!  adhoc_bond(:) = .false.

  return
end subroutine find_neigh

