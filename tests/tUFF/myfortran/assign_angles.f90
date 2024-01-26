subroutine assign_angles
  use conf, only: natoms, neighbors, connectivity, nangles, angles, ufftypes, tol, nangletypes, angletypesmap, bond_orders
  use uff, only: theta_0, Zstar_I
  use pars, only: angle_pars, angle_tags
  implicit none
  integer :: i, ia, ib1, ib2, i1, i2, i3, j, j1, j2, k, m
  real(8) :: r_IJ, r_JK, r_IK, c0, s0
  character(5) :: ati, atj, atk
  logical :: found
  character(5), allocatable :: tmp(:,:)
  real(8), allocatable :: tmp2(:,:)
  
  ! compute the total number of angles
  nangles = 0
  do i2 = 1, natoms
     do j1 = 1, neighbors(i2)-1
        do j2 = j1+1, neighbors(i2)
           nangles = nangles + 1
        end do
     end do
  end do
  
  ! populate angles array
  allocate ( angles(nangles,3), stat = i )
  if ( i /= 0 ) stop 'Allocation failed for ANGLES'
  nangles = 0
  do i2 = 1, natoms
     do j1 = 1, neighbors(i2)-1
        i1 = connectivity(i2,j1)
        do j2 = j1+1, neighbors(i2)
           i3 = connectivity(i2,j2)
           nangles = nangles + 1
           angles(nangles,1) = i1
           angles(nangles,2) = i2
           angles(nangles,3) = i3
        end do
     end do
  end do

  ! calculate number of different angle types and store type arrays
  allocate ( angletypesmap(nangles), stat = i )
  if ( i /= 0 ) stop 'Allocation failed for ANGLETYPESMAP'
  allocate ( tmp(nangles,3), stat = i )
  if ( i /= 0 ) stop 'allocation TMP'
  allocate ( tmp2(nangles,2), stat = i )
  if ( i /= 0 ) stop 'allocation TMP2'
  nangletypes = 0
  do ia = 1, nangles
     call find_bond ( angles(ia,1), angles(ia,2), ib1 )
     call find_bond ( angles(ia,2), angles(ia,3), ib2 )
!write(0,*) ia, ' atom1 ', ufftypes(angles(ia,1)), angles(ia,1)
!write(0,*) ia, ' atom2 ', ufftypes(angles(ia,2)), angles(ia,2)
!write(0,*) ia, ' atom3 ', ufftypes(angles(ia,3)), angles(ia,3)
!write(0,*) ia, ' bond1 ', ib1, bond_orders(ib1)
!write(0,*) ia, ' bond2 ', ib2, bond_orders(ib2)
     found = .false.
     do m = 1, nangletypes
        if ( ufftypes(angles(ia,2)) == tmp(m,2) .and. &
             ( ( ufftypes(angles(ia,1)) == tmp(m,1) .and. abs(bond_orders(ib1)-tmp2(m,1)) < tol .and. &
             ufftypes(angles(ia,3)) == tmp(m,3) .and. abs(bond_orders(ib2)-tmp2(m,2)) < tol ) .or. &
             ( ufftypes(angles(ia,1)) == tmp(m,3) .and. abs(bond_orders(ib1)-tmp2(m,2)) < tol .and. &
             ufftypes(angles(ia,3)) == tmp(m,1) .and. abs(bond_orders(ib2)-tmp2(m,1)) < tol ) ) ) then
           found = .true.
           angletypesmap(ia) = m
!if(ia==7)write(0,*) 'found ', m           
!if(ia==7)write(0,*) 'my ', ufftypes(angles(ia,1)), ufftypes(angles(ia,2)), ufftypes(angles(ia,3))
!if(ia==7)write(0,*) 'old ', tmp(m,1), tmp(m,2), tmp(m,3)
!if(ia==7)write(0,*) 'if ', ufftypes(angles(ia,2)), tmp(m,2)
!if(ia==7)write(0,*) 
!if(ia==7)write(0,*) 
           exit
        end if
     end do
     if ( .not. found ) then
        nangletypes = nangletypes + 1
        tmp(nangletypes,1:3) = ufftypes(angles(ia,1:3))
        tmp2(nangletypes,1) = bond_orders(ib1)
        tmp2(nangletypes,2) = bond_orders(ib2)
        angletypesmap(ia) = nangletypes
     end if
!write(0,*) ia, ' type ', angletypesmap(ia)
!write(0,*)     
  end do
!stop
  
  ! compute pars
  allocate ( angle_pars(nangletypes,4), stat = i )
  if ( i /= 0 ) stop 'Allocation failed for ANGLE_PARS'
  allocate ( angle_tags(nangletypes), stat = i )
  if ( i /= 0 ) stop 'Allocation failed for ANGLE_TAGS'
  do m = 1, nangletypes
     ati = tmp(m,1)
     atj = tmp(m,2)
     atk = tmp(m,3)
     call find_I ( ati, i )
     call find_I ( atj, j )
     call find_I ( atk, k )
     call calc_r_IJ ( i, j, tmp2(m,1), r_IJ )
     call calc_r_IJ ( j, k, tmp2(m,2), r_JK )
     r_IK = sqrt( r_IJ*r_IJ + r_JK*r_JK - 2.d0 * r_IJ * r_JK * cos(theta_0(j)) )
     angle_pars(m,1) = 664.12d0 * Zstar_I(i) * Zstar_I(k) / r_IK**5 * ( 3.d0 * r_IJ * r_JK * sin(theta_0(j))*sin(theta_0(j)) - r_IK*r_IK * cos(theta_0(j)) )
!r_IK = sqrt( r_IJ*r_IJ + r_JK*r_JK - 2.d0 * r_IJ * r_JK * cos(120.d0/180.d0*acos(-1.d0)) )
!angle_pars(m,1) = 664.12d0 * Zstar_I(i) * Zstar_I(k) / r_IK**5 * ( 3.d0 * r_IJ * r_JK * sin(120.d0/180.d0*acos(-1.d0)) - r_IK*r_IK * cos(120.d0/180.d0*acos(-1.d0)) )
     if ( atj(3:3) == '1' .or. atj(3:3) == '2' .or. atj(3:3) == 'R' .or. atj(3:3) == '4' .or. atj(3:3) == '6' ) then ! cosine/periodic
        angle_tags(m) = 'cosine/periodic'
        angle_pars(m,1) = 0.5d0 * angle_pars(m,1)
        angle_pars(m,4) = 0.d0
        if ( atj(3:3) == '1' ) then
           angle_pars(m,2) = 1.d0
           angle_pars(m,3) = 1.d0
        else if ( atj(3:3) == '2' .or. atj(3:3) == 'R' ) then
           angle_pars(m,2) = -1.d0
           angle_pars(m,3) = 3.d0
        else if ( atj(3:3) == '4' .or. atj(3:3) == '6' ) then
           angle_pars(m,2) = 1.d0
           angle_pars(m,3) = 4.d0
        end if
     else if ( atj(3:3) == '3' .or. atj(3:3) == '5' ) then ! fourier
        angle_tags(m) = 'fourier'
        s0 = sin(theta_0(j))
        c0 = cos(theta_0(j))
        angle_pars(m,4) = 1.d0 / ( 4.d0 * s0*s0 )
        angle_pars(m,3) = -4.d0 * angle_pars(m,4) * c0
        angle_pars(m,2) = angle_pars(m,4) * ( 2.d0 * c0*c0 + 1.d0 )
     else
        write(0,*) 'ERROR ASSIGNANGLES: coordination environment not found for atoms ', ati, atj, atk
        write(0,*) 'STOP'
        stop
     end if
     write(*,'(a,a,1x,i4,3(1x,a))') 'angle type ', angle_tags(m), m, trim(ati), trim(atj), trim(atk)
  end do
  deallocate ( tmp, stat = i )
  if ( i /= 0 ) stop 'deallocation TMP'
  deallocate ( tmp2, stat = i )
  if ( i /= 0 ) stop 'deallocation TMP2'
  
  return
end subroutine assign_angles
