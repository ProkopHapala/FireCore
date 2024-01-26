subroutine find_rings
  use conf, only: natoms, neighbors, connectivity
  implicit none
  integer :: ia, nb(0:7), jn1, jn2, jn3, jn4, jn5, jn6
  logical :: explored(natoms), changed, aromatic

  explored(:) = .false.

  do
     changed = .false.
  
     do ia = 1, natoms 
        nb(1) = ia ! 1st atom
        do jn1 = 1, neighbors(nb(1))
           nb(2) = connectivity(nb(1),jn1) ! 2nd atom
           do jn2 = 1, neighbors(nb(2))
              nb(3) = connectivity(nb(2),jn2) ! 3rd atom
              if ( nb(3) == nb(1) ) cycle
              do jn3 = 1, neighbors(nb(3))
                 nb(4) = connectivity(nb(3),jn3) ! 4th atom
                 if ( nb(4) == nb(2) ) cycle
                 do jn4 = 1, neighbors(nb(4))
                    nb(5) = connectivity(nb(4),jn4) ! 5th atom
                    if ( nb(5) == nb(3) ) cycle
                    do jn5 = 1, neighbors(nb(5))
                       nb(6) = connectivity(nb(5),jn5) ! 6th atom
                       if ( nb(6) == nb(4) ) cycle
                       
                       ! found 5-member ring
                       if ( nb(6) == nb(1) ) then
                          if ( explored(nb(1)) .and. explored(nb(2)) .and. explored(nb(3)) .and. explored(nb(4)) .and. explored(nb(5)) ) cycle
                          nb(0) = nb(5)
                          ! count double bonds & heteroatoms
                          call check_aroma ( 5, nb(0:6), aromatic )
                          if ( aromatic ) then
                             changed = .true.
                             explored(nb(1:5)) = .true.
                             call set_aroma ( 5, nb(1:6) )  
                          end if
                       end if

                       do jn6 = 1, neighbors(nb(6))
                          nb(7) = connectivity(nb(6),jn6) ! 7th atom
                          if ( nb(7) == nb(5) ) cycle

                          ! found 6-member ring
                          if ( nb(7) == nb(1) ) then
                             if ( explored(nb(1)) .and. explored(nb(2)) .and. explored(nb(3)) .and. explored(nb(4)) .and. explored(nb(5)) .and. explored(nb(6)) ) cycle
                             nb(0) = nb(6)
                             ! count double bonds & heteroatoms
                             call check_aroma ( 6, nb(0:7), aromatic )
                             if ( aromatic ) then
                                changed = .true.
                                explored(nb(1:6)) = .true.
                                call set_aroma ( 6, nb(1:7) )
                             end if
                          end if
                          
                       end do ! jn6
                    end do ! jn5
                 end do ! jn4
              end do ! jn3
           end do ! jn2
        end do ! jn1
     end do ! ia
     
     if ( .not. changed ) exit
  end do

  return
end subroutine find_rings
