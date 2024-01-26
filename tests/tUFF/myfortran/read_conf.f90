subroutine read_conf ( input )
  use conf, only: natoms, side, side_lammps, tipo, pos, pos_lammps, q, &
       ufftypes, valence_target, neighbors, connectivity, set_atom!, adhoc_atom
  implicit none
  character(200), intent(in) :: input
  real(8), parameter :: tiny = 1.d-10
  integer :: i, j, k
  real(8) :: tmp(3), A(3), B(3), C(3), An(3), pos_tmp(3)
  real(8) :: vec(3), R(3,3), RR(3,3), vol
  character(3) :: comment
  
  ! read xyz
  open(unit=10,file=input)
  read(10,*) natoms
  allocate ( tipo(natoms), stat = i )
  if ( i /= 0 ) stop 'allocation TIPO'
  allocate ( pos(natoms,3), stat = i )
  if ( i /= 0 ) stop 'allocation POS'
  allocate ( pos_lammps(natoms,3), stat = i )
  if ( i /= 0 ) stop 'allocation POS_LAMMPS'
  allocate ( q(natoms), stat = i )
  if ( i /= 0 ) stop 'allocation Q'
  read(10,*) comment, ((side(i,j),j=1,3),i=1,3)
  do i = 1, natoms
     read(10,*) tipo(i), (pos(i,j),j=1,3), q(i)
     if ( tipo(i) == 'h' ) tipo(i) = 'H'
     if ( tipo(i) == 'c' ) tipo(i) = 'C'
     if ( tipo(i) == 'n' ) tipo(i) = 'N'
     if ( tipo(i) == 'o' ) tipo(i) = 'O'
  end do
  close(10)

  side_lammps(:,:) = side(:,:)
  pos_lammps(:,:) = pos(:,:)
  
  ! check if the cell vectors form a right-handed basis
  call cross_product ( side(1,:), side(2,:), vec )
  if ( dot_product( vec, side(3,:) ) < 0.d0 ) then
     write(*,*) 'WARNING READCONF: cell vectors are left-handed!'
     write(*,*) 'WARNING READCONF: b and c cell vectors will be swapped.'
     do i = 1, 3
        side_lammps(3,i) = side(2,i)
        side_lammps(2,i) = side(3,i)
     end do
     do i = 1, natoms
        pos_lammps(i,3) = pos(i,2)
        pos_lammps(i,2) = pos(i,3)
     end do
  end if
  
  ! see https://docs.lammps.org/Howto_triclinic.html
  A(:) = side_lammps(1,:)
  tmp(1) = sqrt( dot_product(A,A) )
  An(:) = A(:) / tmp(1)
  B(:) = side_lammps(2,:)
  C(:) = side_lammps(3,:)
  side_lammps(1,1) = sqrt( dot_product(A,A) )
  side_lammps(1,2:3) = 0.d0
  side_lammps(2,1) = dot_product(B,An)
  side_lammps(2,2) = sqrt( dot_product(B,B) - side_lammps(2,1)*side_lammps(2,1) )
  side_lammps(2,3) = 0.d0
  side_lammps(3,1) = dot_product(C,An)
  side_lammps(3,2) = ( dot_product(B,C) - side_lammps(2,1)*side_lammps(3,1) ) / side_lammps(2,2)
  side_lammps(3,3) = sqrt( dot_product(C,C) - side_lammps(3,1)*side_lammps(3,1) - side_lammps(3,2)*side_lammps(3,2) )
  call cross_product ( B, C, vec )
  vol = dot_product( A, vec )
  do i = 1, 3
     R(1,i) = vec(i) / vol
  end do
  call cross_product ( C, A, vec )
  do i = 1, 3
     R(2,i) = vec(i) / vol
  end do
  call cross_product ( A, B, vec )
  do i = 1, 3
     R(3,i) = vec(i) / vol
  end do
  RR = matmul( side_lammps, R )

  ! final check for lammps settings
  if ( abs(side_lammps(1,2)) > tiny .or. abs(side_lammps(1,3)) > tiny .or. abs(side_lammps(2,3)) > tiny ) then
     write(0,*) 'ERROR READCONF something went wrong...'
     write(0,*) 'ERROR READCONF cell after rotation:'
     write(0,*) side_lammps(1,:)
     write(0,*) side_lammps(2,:)
     write(0,*) side_lammps(3,:)
     write(0,*) 'STOP'
     stop
  end if

  ! rotate atoms according to the new cell orientation
  do k = 1, natoms
     pos_tmp(:) = pos_lammps(k,:)
     pos_lammps(k,:) = 0.d0
     do i = 1, 3
        do j = 1, 3
           pos_lammps(k,i) = pos_lammps(k,i) + RR(i,j) * pos_tmp(j)
        end do
     end do
  end do

  ! allocate and initialize atomic arrays
  allocate ( ufftypes(natoms), stat = i )
  if ( i /= 0 ) stop 'allocation UFFTYPES'
  ufftypes(:) = ''
  allocate ( valence_target(natoms), stat = i )
  if ( i /= 0 ) stop 'allocation VALENCE_TARGET'
  do i = 1, natoms
     if ( tipo(i) == 'H' ) valence_target(i) = 1
     if ( tipo(i) == 'C' ) valence_target(i) = 4
     if ( tipo(i) == 'N' ) valence_target(i) = 3
     if ( tipo(i) == 'O' ) valence_target(i) = 2
  end do
  allocate ( neighbors(natoms), stat = i )
  if ( i /= 0 ) stop 'allocation NEIGHBORS'
  neighbors(:) = 0
  allocate ( connectivity(natoms,4), stat = i )
  if ( i /= 0 ) stop 'allocation CONNECTIVITY'
  connectivity(:,:) = -1
  allocate ( set_atom(natoms), stat = i )
  if ( i /= 0 ) stop 'allocation SET_ATOM'
  set_atom(:) = .false.
!  allocate ( adhoc_atom(natoms), stat = i )
!  if ( i /= 0 ) stop 'allocation ADHOC_ATOM'
!  adhoc_atom(:) = .false.
  
  return
end subroutine read_conf
