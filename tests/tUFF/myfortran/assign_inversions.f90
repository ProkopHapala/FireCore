subroutine assign_inversions
  use conf, only: natoms, neighbors, connectivity, ninversions, inversions, ufftypes, ninversiontypes, inversiontypesmap
  use pars, only: inversion_pars
  implicit none
  real(8), parameter :: pi = acos(-1.d0)
  integer :: i, ii, i1, i2, i3, i4, m
  real(8) :: omega0 = 0.d0
  character(5) :: ati, atj, atk, atl
  logical :: found
  character(5), allocatable :: tmp(:,:)
  
  ! compute the total number of improper torsions
  ninversions = 0
  do i1 = 1, natoms
     if ( neighbors(i1) == 3 ) then
        ninversions = ninversions + 3
     end if
  end do
  
  ! populate improper torsions array
  allocate ( inversions(ninversions,4), stat = i )
  if ( i /= 0 ) stop 'Allocation failed for INVERSIONS'
  ninversions = 0
  do i1 = 1, natoms
     if ( neighbors(i1) == 3 ) then
        i2 = connectivity(i1,1)
        i3 = connectivity(i1,2)
        i4 = connectivity(i1,3)
        ninversions = ninversions + 1
        inversions(ninversions,1) = i1
        inversions(ninversions,2) = i2
        inversions(ninversions,3) = i3
        inversions(ninversions,4) = i4
        ninversions = ninversions + 1
        inversions(ninversions,1) = i1
        inversions(ninversions,2) = i2
        inversions(ninversions,3) = i4
        inversions(ninversions,4) = i3
        ninversions = ninversions + 1
        inversions(ninversions,1) = i1
        inversions(ninversions,2) = i4
        inversions(ninversions,3) = i3
        inversions(ninversions,4) = i2
     end if
  end do

  ! calculate number of different improper torsion types and store type arrays
  allocate ( inversiontypesmap(ninversions), stat = i )
  if ( i /= 0 ) stop 'Allocation failed for INVERSIONTYPESMAP'
  allocate ( tmp(ninversions,4), stat = i )
  if ( i /= 0 ) stop 'allocation TMP'
  ninversiontypes = 0
  do ii = 1, ninversions
     found = .false.
     do m = 1, ninversiontypes
        if ( ufftypes(inversions(ii,1)) == tmp(m,1) .and. &
             ( ufftypes(inversions(ii,2)) == tmp(m,2) .and. ufftypes(inversions(ii,3)) == tmp(m,3) .and. ufftypes(inversions(ii,4)) == tmp(m,4) ) .or. &
             ( ufftypes(inversions(ii,2)) == tmp(m,2) .and. ufftypes(inversions(ii,3)) == tmp(m,4) .and. ufftypes(inversions(ii,4)) == tmp(m,3) ) .or. &
             ( ufftypes(inversions(ii,2)) == tmp(m,3) .and. ufftypes(inversions(ii,3)) == tmp(m,2) .and. ufftypes(inversions(ii,4)) == tmp(m,4) ) .or. &
             ( ufftypes(inversions(ii,2)) == tmp(m,3) .and. ufftypes(inversions(ii,3)) == tmp(m,4) .and. ufftypes(inversions(ii,4)) == tmp(m,2) ) .or. &
             ( ufftypes(inversions(ii,2)) == tmp(m,4) .and. ufftypes(inversions(ii,3)) == tmp(m,2) .and. ufftypes(inversions(ii,4)) == tmp(m,3) ) .or. &
             ( ufftypes(inversions(ii,2)) == tmp(m,4) .and. ufftypes(inversions(ii,3)) == tmp(m,3) .and. ufftypes(inversions(ii,4)) == tmp(m,2) ) ) then
           found = .true.
           inversiontypesmap(ii) = m
           exit
        end if
     end do
     if ( .not. found ) then
        ninversiontypes = ninversiontypes + 1
        tmp(ninversiontypes,1:4) = ufftypes(inversions(ii,1:4))
        inversiontypesmap(ii) = ninversiontypes
     end if
  end do
  
  ! compute pars
  allocate ( inversion_pars(ninversiontypes,4), stat = i )
  if ( i /= 0 ) stop 'Allocation failed for INVERSION_PARS'
  do m = 1, ninversiontypes
     ati = tmp(m,1)
     atj = tmp(m,2)
     atk = tmp(m,3)
     atl = tmp(m,4)
     ! sp2 carbon
     if ( ati(1:3) == 'C_2' .or. ati(1:3) == 'C_R' ) then 
        ! carbonyl
        if ( atj(1:3) == 'O_2' .or. atk(1:3) == 'O_2' .or. atl(1:3) == 'O_2' ) then
           inversion_pars(m,1) = 50.d0  
        else
           inversion_pars(m,1) = 6.d0
        end if
        inversion_pars(m,2) = 1.d0
        inversion_pars(m,3) = -1.d0
        inversion_pars(m,4) = 0.d0
     ! sp2 nitrogen
     else if ( ati(1:3) == 'N_2' .or. ati(1:3) == 'N_R' ) then
        inversion_pars(m,1) = 6.d0
        inversion_pars(m,2) = 1.d0
        inversion_pars(m,3) = -1.d0
        inversion_pars(m,4) = 0.d0
     ! sp3 nitrogen
     else if ( ati(1:3) == 'N_3' ) then
        inversion_pars(m,1) = 0.d0
        inversion_pars(m,2) = 0.d0
        inversion_pars(m,3) = 0.d0
        inversion_pars(m,4) = 0.d0
     ! sp3 group 15
     else if ( ati(1:3) == 'P_3' .or. ati(1:3) == 'As3' .or. ati(1:3) == 'Sb3' .or. ati(1:3) == 'Bi3' ) then
        if ( ati(1:3) == 'P_3' ) then
           omega0 = 84.4339d0 / 180.d0 * pi
        else if ( ati(1:3) == 'As3' ) then
           omega0 = 86.9735d0 / 180.d0 * pi
        else if ( ati(1:3) == 'Sb3' ) then
           omega0 = 87.7047d0 / 180.d0 * pi
        else if ( ati(1:3) == 'Bi3' ) then
           omega0 = 90.d0 / 180.d0 * pi
        end if
        inversion_pars(m,2) = 4.d0 * cos(omega0)*cos(omega0) - cos(2.d0*omega0)
        inversion_pars(m,3) = -4.d0 * cos(omega0)
        inversion_pars(m,4) = 1.d0
        inversion_pars(m,1) = 22.d0 / ( inversion_pars(m,2) + inversion_pars(m,3) + inversion_pars(m,4) )
     else
        write(0,*) 'ERROR ASSIGNINVERSIONS: inversion case not found for atoms ', ati, atj, atk, atl
        write(0,*) 'STOP'
        stop
     end if
     inversion_pars(m,1) = inversion_pars(m,1) / 3.d0
     write(*,'(a,i4,4(1x,a))') 'inversion type ', m, trim(ati), trim(atj), trim(atk), trim(atl)
  end do
  deallocate ( tmp, stat = i )
  if ( i /= 0 ) stop 'deallocation TMP'

  return
end subroutine assign_inversions
