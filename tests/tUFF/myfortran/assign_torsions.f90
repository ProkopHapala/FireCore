subroutine assign_torsions
  use conf, only: natoms, neighbors, connectivity, ntorsions, torsions, ufftypes, tol, ntorsiontypes, torsiontypesmap, bond_orders
  use uff, only: V_I, U_I
  use pars, only: torsion_pars
  implicit none
  integer :: i, ib, it, i1, i2, i3, i4, j, j1, j2, j3, k, m!, ntors
  character(5) :: ati, atj, atk, atl
  logical :: found
  character(5), allocatable :: tmp(:,:)
  real(8), allocatable :: tmp2(:)
  integer, allocatable :: tmp3(:), tmp4(:,:), tmp5(:,:)
  
  ! compute the total number of proper torsions
  ntorsions = 0
  do i1 = 1, natoms
     do j1 = 1, neighbors(i1)
        i2 = connectivity(i1,j1)
        do j2 = 1, neighbors(i2)
           i3 = connectivity(i2,j2)   
           if ( i3 /= i1 ) then
              do j3 = 1, neighbors(i3)
                 i4 = connectivity(i3,j3)   
                 if ( i4 /= i2 .and. i4 > i1 ) then ! avoid 3-membered rings and double-counting
                    if ( ufftypes(i2)(3:3) == '1' .or. ufftypes(i3)(3:3) == '1' ) cycle ! Torsional potentials for central bonds involving sp-hybridized centers X-1 were assigned a value of zero
                    ntorsions = ntorsions + 1
                 end if
              end do
           end if
        end do
     end do
  end do
  
  ! populate proper torsions array
  allocate ( torsions(ntorsions,4), stat = i )
  if ( i /= 0 ) stop 'Allocation failed for TORSIONS'
  ntorsions = 0
  do i1 = 1, natoms
     do j1 = 1, neighbors(i1)
        i2 = connectivity(i1,j1)
        do j2 = 1, neighbors(i2)
           i3 = connectivity(i2,j2)   
           if ( i3 /= i1 ) then
              do j3 = 1, neighbors(i3)
                 i4 = connectivity(i3,j3)   
                 if ( i4 /= i2 .and. i4 > i1 ) then ! avoid 3-membered rings and double-counting
                    if ( ufftypes(i2)(3:3) == '1' .or. ufftypes(i3)(3:3) == '1' ) cycle ! Torsional potentials for central bonds involving sp-hybridized centers X-1 were assigned a value of zero
                    ntorsions = ntorsions + 1
                    torsions(ntorsions,1) = i1
                    torsions(ntorsions,2) = i2
                    torsions(ntorsions,3) = i3
                    torsions(ntorsions,4) = i4
                 end if
              end do
           end if
        end do
     end do
  end do

  ! calculate number of different torsion types and store type arrays
  allocate ( torsiontypesmap(ntorsions), stat = i )
  if ( i /= 0 ) stop 'Allocation failed for TORSIONTYPESMAP'
  allocate ( tmp(ntorsions,4), stat = i )
  if ( i /= 0 ) stop 'allocation TMP'
  allocate ( tmp2(ntorsions), stat = i )
  if ( i /= 0 ) stop 'allocation TMP2'
  allocate ( tmp3(ntorsions), stat = i )
  if ( i /= 0 ) stop 'Allocation failed for TMP3'
  allocate ( tmp4(ntorsions,2), stat = i )
  if ( i /= 0 ) stop 'Allocation failed for TMP4'
  allocate ( tmp5(ntorsions,4), stat = i )
  if ( i /= 0 ) stop 'Allocation failed for TMP5'
  ntorsiontypes = 0
  do it = 1, ntorsions
     call find_bond ( torsions(it,2), torsions(it,3), ib )
     found = .false.
     do m = 1, ntorsiontypes
!        if ( abs(bond_orders(ib)-tmp2(m)) < tol .and. &
!             ( ufftypes(torsions(it,1)) == tmp(m,1) .and. ufftypes(torsions(it,2)) == tmp(m,2) .and. &
!             ufftypes(torsions(it,3)) == tmp(m,3) .and. ufftypes(torsions(it,4)) == tmp(m,4) ) .or. &
!             ( ufftypes(torsions(it,1)) == tmp(m,4) .and. ufftypes(torsions(it,2)) == tmp(m,3) .and. &
!             ufftypes(torsions(it,3)) == tmp(m,2) .and. ufftypes(torsions(it,4)) == tmp(m,1) ) ) then
        if ( abs(bond_orders(ib)-tmp2(m)) > tol ) cycle ! wrong bond order
        if ( ( ufftypes(torsions(it,1)) == tmp(m,1) .and. ufftypes(torsions(it,2)) == tmp(m,2) .and. ufftypes(torsions(it,3)) == tmp(m,3) .and. ufftypes(torsions(it,4)) == tmp(m,4) ) .or. &
             ( ufftypes(torsions(it,1)) == tmp(m,4) .and. ufftypes(torsions(it,2)) == tmp(m,3) .and. ufftypes(torsions(it,3)) == tmp(m,2) .and. ufftypes(torsions(it,4)) == tmp(m,1) ) ) then
           if ( ( ufftypes(torsions(it,1)) == tmp(m,1) .and. neighbors(torsions(it,2)) == tmp4(m,1) .and. neighbors(torsions(it,3)) == tmp4(m,2) ) .or. &
                ( ufftypes(torsions(it,1)) == tmp(m,4) .and. neighbors(torsions(it,2)) == tmp4(m,2) .and. neighbors(torsions(it,3)) == tmp4(m,1) ) ) then ! check neighbors
              found = .true.
              torsiontypesmap(it) = m
              exit
           end if
        end if
     end do
     if ( .not. found ) then
        ntorsiontypes = ntorsiontypes + 1
        tmp(ntorsiontypes,1:4) = ufftypes(torsions(it,1:4))
        tmp2(ntorsiontypes) = bond_orders(ib)
        torsiontypesmap(it) = ntorsiontypes
        tmp3(ntorsiontypes) = ( neighbors(torsions(it,2)) - 1 ) * ( neighbors(torsions(it,3)) - 1 )
        tmp4(ntorsiontypes,1) = neighbors(torsions(it,2))
        tmp4(ntorsiontypes,2) = neighbors(torsions(it,3))
        tmp5(ntorsiontypes,1:4) = torsions(it,1:4)
     end if
  end do

  ! compute pars
  allocate ( torsion_pars(ntorsiontypes,3), stat = i )
  if ( i /= 0 ) stop 'Allocation failed for TORSION_PARS'
  do m = 1, ntorsiontypes
     ati = tmp(m,1)
     atj = tmp(m,2)
     atk = tmp(m,3)
     atl = tmp(m,4)
     call find_I ( atj, j )
     call find_I ( atk, k )

     ! number of torsions
!     ntors = 1
!     if ( atj(1:2) == 'N_' .or. atj(1:2) == 'P_' .or. atj(1:2) == 'As' .or. atj(1:2) == 'Sb' .or. atj(1:2) == 'Bi' ) then
!        if ( atj(3:3) == '3' ) ntors = ntors * 2 ! e.g. sp3 N -> 2 neighs
!     else if ( atj(1:2) == 'C_' .or. atj(1:2) == 'Si' .or. atj(1:2) == 'Ge' .or. atj(1:2) == 'Sn' .or. atj(1:2) == 'Pb' ) then
!        if ( atj(3:3) == '3' ) then
!           ntors = ntors * 3 ! e.g. sp3 C -> 3 neighs
!        else if ( atj(3:3) == '2' .or. atj(3:3) == 'R' ) then
!           ntors = ntors * 2 ! e.g. sp2 C -> 2 neighs
!        end if
!     end if
!     if ( atk(1:2) == 'N_' .or. atk(1:2) == 'P_' .or. atk(1:2) == 'As' .or. atk(1:2) == 'Sb' .or. atk(1:2) == 'Bi' ) then
!        if ( atk(3:3) == '3' ) ntors = ntors * 2 ! e.g. sp3 N -> 2 neighs
!     else if ( atk(1:2) == 'C_' .or. atk(1:2) == 'Si' .or. atk(1:2) == 'Ge' .or. atk(1:2) == 'Sn' .or. atk(1:2) == 'Pb' ) then
!        if ( atk(3:3) == '3' ) then
!           ntors = ntors * 3 ! e.g. sp3 C -> 3 neighs
!        else if ( atk(3:3) == '2' .or. atk(3:3) == 'R' ) then
!           ntors = ntors * 2 ! e.g. sp2 C -> 2 neighs
!        end if
!     end if
     
     ! case (a) * - sp3 - sp3 - *
     if ( atj(3:3) == '3' .and. atk(3:3) == '3' ) then 
        ! specific general case
        torsion_pars(m,1) = sqrt( V_I(j) * V_I(k) )
        torsion_pars(m,2) = 1.d0
        torsion_pars(m,3) = 3.d0
        ! special case of * - group 16 sp3 - group 16 sp3 - *
        if ( ( atj(1:2) == 'O_' .or. atj(1:2) == 'S_' .or. atj(1:2) == 'Se' .or. atj(1:2) == 'Te' .or. atj(1:2) == 'Po' ) .and. &
             ( atk(1:2) == 'O_' .or. atk(1:2) == 'S_' .or. atk(1:2) == 'Se' .or. atk(1:2) == 'Te' .or. atk(1:2) == 'Po' ) ) then 
           torsion_pars(m,3) = 2.d0
           torsion_pars(m,1) = 1.d0
           if ( atj(1:2) == 'O_' ) then
              torsion_pars(m,1) = torsion_pars(m,1) * 2.d0
           else
              torsion_pars(m,1) = torsion_pars(m,1) * 6.8d0
           end if
           if ( atk(1:2) == 'O_' ) then
              torsion_pars(m,1) = torsion_pars(m,1) * 2.d0
           else
              torsion_pars(m,1) = torsion_pars(m,1) * 6.8d0
           end if
           torsion_pars(m,1) = sqrt( torsion_pars(m,1) )
        end if
     
     ! (b) * - sp3 - sp2 - *
     else if ( ( atj(3:3) == '3' .and. ( atk(3:3) == '2' .or. atk(3:3) == 'R' ) ) .or. &
          ( ( atj(3:3) == '2' .or. atj(3:3) == 'R' ) .and. atk(3:3) == '3' ) ) then 
        ! specific general case
        torsion_pars(m,1) = 1.d0
        torsion_pars(m,2) = -1.d0
        torsion_pars(m,3) = 6.d0
        ! special case of * - group 16 sp3 - sp2 - *
        if ( ( atj(3:3) == '3' .and. ( atj(1:2) == 'O_' .or. atj(1:2) == 'S_' .or. atj(1:2) == 'Se' .or. atj(1:2) == 'Te' .or. atj(1:2) == 'Po' ) ) .or. &
             ( atk(3:3) == '3' .and. ( atk(1:2) == 'O_' .or. atk(1:2) == 'S_' .or. atk(1:2) == 'Se' .or. atk(1:2) == 'Te' .or. atk(1:2) == 'Po' ) ) ) then
           torsion_pars(m,1) = 5.d0 * sqrt( U_I(j) * U_I(k) ) * ( 1.d0 + 4.18d0 * log( tmp2(m) ) )
           torsion_pars(m,2) = 1.d0
           torsion_pars(m,3) = 2.d0
!        ! special case of * - sp3 - sp2 - sp2
!        else if ( ( atj(3:3) == '3' .and. ( atl(3:3) == '2' .or. atl(3:3) == 'R' ) ) .or. ( atk(3:3) == '3' .and. ( ati(3:3) == '2' .or. ati(3:3) == 'R' ) ) ) then
!           torsion_pars(m,1) = 2.d0
!           torsion_pars(m,2) = 1.d0
!           torsion_pars(m,3) = 3.d0
        end if

        ! special case of * - sp3 - sp2 bounded to another sp2 atom
        if ( atj(3:3) == '3' ) then
           found = .false.
           do i = 1, neighbors(tmp5(m,3))
              j = connectivity(tmp5(m,3),i)
              if ( ufftypes(j)(3:3) == '2' .or. ufftypes(j)(3:3) == 'R' ) then
                 found = .true.
                 exit
              end if
           end do
           if ( found ) then
              torsion_pars(m,1) = 2.d0
              torsion_pars(m,2) = 1.d0
              torsion_pars(m,3) = 3.d0
           end if
        else if ( atk(3:3) == '3' ) then
           found = .false.
           do i = 1, neighbors(tmp5(m,2))
              j = connectivity(tmp5(m,2),i)
              if ( ufftypes(j)(3:3) == '2' .or. ufftypes(j)(3:3) == 'R' ) then
                 found = .true.
                 exit
              end if
           end do
           if ( found ) then
              torsion_pars(m,1) = 2.d0
              torsion_pars(m,2) = 1.d0
              torsion_pars(m,3) = 3.d0
           end if
        end if
        
     ! (c) * - sp2 - sp2 - *
     else if ( ( atj(3:3) == '2' .or. atj(3:3) == 'R' ) .and. ( atk(3:3) == '2' .or. atk(3:3) == 'R' ) ) then 
        ! specific general case
        torsion_pars(m,1) = 5.d0 * sqrt( U_I(J) * U_I(K) ) * ( 1.d0 + 4.18d0 * log( tmp2(m) ) )
        torsion_pars(m,2) = -1.d0
        torsion_pars(m,3) = 2.d0
        ! special case for carboxyl from Prokop
!        if ( ( ati == 'H_' .and. atj == 'O_R' .and. atk == 'C_2' ) .or. ( atl == 'H_' .and. atk == 'O_R' .and. atj == 'C_2' ) ) torsion_pars(m,3) = 1.d0

!     ! special case of * - sp - * - *
!     else if ( atj(3:3) == '1' .or. atk(3:3) == '1' ) then
!        torsion_pars(m,1) = 0.d0
!        torsion_pars(m,2) = 1.d0
!        torsion_pars(m,3) = 1.d0

     else
        write(0,*) 'ERROR ASSIGNTORSIONS: torsion case not found for atoms ', ati, atj, atk, atl
        write(0,*) 'STOP'
        stop
     end if
!     torsion_pars(m,1) = 0.5d0 * torsion_pars(m,1) / dble(ntors)
     torsion_pars(m,1) = 0.5d0 * torsion_pars(m,1) / dble(tmp3(m))
!if(ati=='H_'.and.atj=='C_R'.and.atk=='N_R'.and.atl=='C_R')write(0,*) 'superscio'     
!if(ati=='C_1'.and.atj=='C_2'.and.atk=='C_3'.and.atl=='H_')then
!   write(0,*) 'superscio'
!   write(0,*) tmp3(m)
!end if
     write(*,'(a,i4,4(1x,a))') 'torsion type ', m, trim(ati), trim(atj), trim(atk), trim(atl)
  end do
  deallocate ( tmp, stat = i )
  if ( i /= 0 ) stop 'deallocation TMP'
  deallocate ( tmp2, stat = i )
  if ( i /= 0 ) stop 'deallocation TMP2'
  deallocate ( tmp3, stat = i )
  if ( i /= 0 ) stop 'deallocation TMP3'
  deallocate ( tmp4, stat = i )
  if ( i /= 0 ) stop 'deallocation TMP4'
  deallocate ( tmp5, stat = i )
  if ( i /= 0 ) stop 'deallocation TMP5'

  return
end subroutine assign_torsions
