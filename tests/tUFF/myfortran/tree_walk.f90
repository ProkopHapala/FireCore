subroutine tree_walk
  use conf, only: natoms, nbonds, bond_orders_int, neighbors, connectivity
  implicit none
  integer, parameter :: nattempt = 3 ! hard coded
  integer :: iattempt, ia1, na1, ib1, ia2, na2, ib2, ia3, na3, ib3
  integer :: bond_orders_tmp(nbonds,0:nattempt)
  logical :: done, exclude
  
  ! see if the problem is trivial
  call check_all ( bond_orders_int, done ) 
  if ( done ) return

  ! --- explore the graph
  ! NB: at this point, only sp2 C, sp2 N and sp C (both triple bond and 2x double bond) are left to be assigned
  ! but triple bonds in the backbone (i.e. cumulene) will be sorted out later

  bond_orders_tmp(:,0) = bond_orders_int(:)
  iattempt = 1
  ! pick an atom to attempt assigning a double bond
  do ia1 = 1, natoms
     call check_atom ( ia1, bond_orders_tmp(:,iattempt-1), exclude )
     if ( exclude ) cycle
     ! pick a a neighbor
     do na1 = 1, neighbors(ia1)
        call find_bond ( ia1, connectivity(ia1,na1), ib1 )
        if ( bond_orders_tmp(ib1,iattempt-1) > 0 ) cycle
        ! try to set it up as a double bond
        bond_orders_tmp(:,iattempt) = bond_orders_tmp(:,iattempt-1)
        bond_orders_tmp(ib1,iattempt) = 2
        ! see if this solves the problem
        call check_all ( bond_orders_tmp(:,iattempt), done ) 
        if ( done ) then
           bond_orders_int(:) = bond_orders_tmp(:,iattempt)
           return
        end if
        
        ! otherwise, continue with the second nested attempt
        iattempt = 2
        ! pick an atom to attempt assigning a double bond
        do ia2 = 1, natoms
           call check_atom ( ia2, bond_orders_tmp(:,iattempt-1), exclude )
           if ( exclude ) cycle
           ! pick a a neighbor
           do na2 = 1, neighbors(ia2)
              call find_bond ( ia2, connectivity(ia2,na2), ib2 )
              if ( bond_orders_tmp(ib2,iattempt-1) > 0 ) cycle
              ! try to set it up as a double bond
              bond_orders_tmp(:,iattempt) = bond_orders_tmp(:,iattempt-1)
              bond_orders_tmp(ib2,iattempt) = 2
              ! see if this solves the problem
              call check_all ( bond_orders_tmp(:,iattempt), done ) 
              if ( done ) then
                 bond_orders_int(:) = bond_orders_tmp(:,iattempt)
                 return
              end if
              
              ! otherwise, continue with the last nested attempt
              iattempt = 3
              ! pick an atom to attempt assigning a double bond
              do ia3 = 1, natoms
                 call check_atom ( ia3, bond_orders_tmp(:,iattempt-1), exclude )
                 if ( exclude ) cycle
                 ! pick a a neighbor
                 do na3 = 1, neighbors(ia3)
                    call find_bond ( ia3, connectivity(ia3,na3), ib3 )
                    if ( bond_orders_tmp(ib3,iattempt-1) > 0 ) cycle
                    ! try to set it up as a double bond
                    bond_orders_tmp(:,iattempt) = bond_orders_tmp(:,iattempt-1)
                    bond_orders_tmp(ib3,iattempt) = 2
                    ! see if this solves the problem
                    call check_all ( bond_orders_tmp(:,iattempt), done ) 
                    if ( done ) then
                       bond_orders_int(:) = bond_orders_tmp(:,iattempt)
                       return
                    end if
                    
                 end do
              end do
              
           end do
        end do
        
     end do
  end do
  
  write(0,*) 'ERROR TREEWALK: tree walk failed'
  write(0,*) 'STOP'
  stop
end subroutine tree_walk
