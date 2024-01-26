subroutine assign_vdw
  use conf, only: natoms, atomtypesmap, natomtypes, ufftypes
  use uff, only: D_I, x_I
  use pars, only: vdw_pars
  implicit none
  real(8), parameter :: fact = 1.d0 / 2.d0**(1.d0/6.d0)
  integer :: i, ia, m
  logical :: found
  character(5) :: tmp(natoms)
  
  ! calculate number of different atom types and store type arrays
  allocate ( atomtypesmap(natoms), stat = i )
  if ( i /= 0 ) stop 'allocation ATOMTYPESMAP'
  natomtypes = 0
  do ia = 1, natoms
     found = .false.
     do m = 1, natomtypes
        if ( ufftypes(ia) == tmp(m) ) then
           found = .true.
           atomtypesmap(ia) = m
           exit
        end if
     end do
     if ( .not. found ) then
        natomtypes = natomtypes + 1
        tmp(natomtypes) = ufftypes(ia)
        atomtypesmap(ia) = natomtypes
     end if
  end do
   
  ! assign vdw parameters
  allocate ( vdw_pars(natomtypes,2), stat = i )
  if ( i /= 0 ) stop 'allocation VDW_PARS'
  do m = 1, natomtypes
     call find_I ( tmp(m), i )
     vdw_pars(m,1) = D_I(i)
     vdw_pars(m,2) = fact * x_I(i)
     write(*,'(a,i4,1x,a)') 'atom type ', m, trim(tmp(m))
  end do

  return
end subroutine assign_vdw
