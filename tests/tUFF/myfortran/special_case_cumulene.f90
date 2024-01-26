subroutine special_case_cumulene
  use conf, only: natoms, bond_orders, set_bond, neighbors, connectivity, &
       bond_orders, bond_orders_int, tipo
  implicit none
  integer :: ia1, ia2, ia3, ia4, ib1, ib2, ib3, ib4, n1, n2, n3

  ! check and manually transform cumulene into triple, double and single bonds
  main: do ia1 = 1, natoms
     if ( tipo(ia1) /= 'C' ) cycle    ! exclude non-C
     if ( neighbors(ia1) /= 2 ) cycle ! looking for an sp C
     do n1 = 1, neighbors(ia1)
        ia2 = connectivity(ia1,n1)
        if ( tipo(ia2) /= 'C' ) cycle    ! exclude non-C
        if ( neighbors(ia2) /= 2 ) cycle ! looking for the second sp C
        call find_bond ( ia1, ia2, ib1 )
        ia4 = connectivity(ia1,mod(n1,2)+1)
        if ( tipo(ia4) /= 'C' ) cycle    ! exclude non-C
        if ( neighbors(ia4) /= 3 ) cycle ! looking for an sp2 C
        call find_bond ( ia1, ia4, ib4 )
        do n2 = 1, neighbors(ia2)
           ia3 = connectivity(ia2,n2)
           if ( ia3 == ia1 ) cycle          ! looking for the other neighbor of ia2
           if ( tipo(ia3) /= 'C' ) cycle    ! exclude non-C
           if ( neighbors(ia3) /= 3 ) cycle ! looking for an sp2 C
           do n3 = 1, neighbors(ia3)
              if ( connectivity(ia3,n3) /= ia4 ) cycle ! found the right ia3
              call find_bond ( ia2, ia3, ib2 )
              call find_bond ( ia3, ia4, ib3 )
              bond_orders(ib1) = 3.d0
              bond_orders(ib2) = 1.d0
              bond_orders(ib3) = 2.d0
              bond_orders(ib4) = 1.d0
              bond_orders_int(ib1) = 3
              bond_orders_int(ib2) = 1
              bond_orders_int(ib3) = 2
              bond_orders_int(ib4) = 1
              set_bond(ib1) = .true.
              set_bond(ib2) = .true.
              set_bond(ib3) = .true.
              set_bond(ib4) = .true.
              cycle main
           end do
        end do
     end do
  end do main

  return
end subroutine special_case_cumulene
