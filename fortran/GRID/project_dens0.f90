
! ===========================================================================
! ===========================================================================
 subroutine project_dens0( ewfaux, f_mul )
   use configuration
   use dimensions
   use interactions
   use neighbor_map
   use grid
   use density
   use charges
   !use outputs
   use options
   implicit none

! Argument Declaration and Description
! ===========================================================================
   !integer iband0,iband1
   !real, dimension (:), pointer, intent (out) :: ewfaux
   real, dimension (nrm), intent (out) :: ewfaux
   ! real, target, dimension (:), allocatable :: ewfaux
   real, intent (in) :: f_mul

! Local Variable Declaration and Description
! ===========================================================================
   integer iatom
   integer imu, inu
   integer in1, in2
   integer mbeta
   integer ineigh
   integer index
   integer index0
   integer ind
   integer i 
   integer issh
   integer imesh
   integer info

   real distX
   real qtot
   real renorm
   real dens

   real, dimension (3) :: r1
   real, dimension (3) :: u
   real, dimension (3) :: dXr
   real, dimension (3) :: g1
   real, dimension (3) :: X0
   real, dimension (3) :: u1X

   real, dimension (numorb_max)   :: psi1

   real, dimension (3,3)          :: lmat
   real, dimension (3,3)          :: invl

!                 + X0 (iatom)
!                / \     u1X = g1 - X0
! uX0 = X0- g0  /   \
!              /     + g1 (nearest grid point to iatom)
!             /
!            /
!           +
!          g0 (origin of grid coords)
!
! g1 ... point where the origin of atomic mesh (sphere) is placed
!  vectors we need:
!    X0g0
!    u1X = g1 - X0
!
! Procedure
! ===========================================================================

!debug
   qtot = 0.0
   renorm = 0.0d0
   lmat = transpose(elvec)   ! copy and invert original elvec to get form written above
   call inv3x3 (lmat,invl)   ! inverse A: solving A*n=x -> n=A-1*x
   do iatom = 1, natoms   ! Loop over atoms
    in1 = imass(iatom)
    r1(:) = ratom(:,iatom)
     u(:) = ratom(:,iatom) - g0(:)
    call mult3x1 (invl,u)   ! get n-vector
    u1X = u
    call mult3x1 (lmat,u1X)
    do i = 1,3
      if (abs(u1X(i)-(ratom(i,iatom) - g0(i))) .gt. 0.000001d0 ) then
         write (*,*) ' ERORR : vec ', u1X(:), 'diff = ', abs(u1X(i)-(ratom(i,iatom) - g0(i)))
         stop
       endif
    enddo
    call mapAtomToGrid( u, X0, g1, index0 )
    u1X(:)           = g1(:) - X0(:)                  ! vector pointing from g1 to X0
    ratom2g(:,iatom) = X0(:)                          ! save iatom coord within the grid unit cell
    ! skip certain atoms to avoid the overcounting; it means all jatoms having the identification number less then iatom because those pairs have been counted previously
    ! Loop over points in the atomic mesh gP
    do imesh = 1, nam
     index = index0 + am2rc(imesh)
      dXr(:) = ram2rc(:,imesh) + u1X(:)
     distX = sqrt(dXr(1)**2 + dXr(2)**2 + dXr(3)**2)  ! distance between the mesh point and iatom
     call getAtomBasis( in1, distX, dXr, psi1 )
     ind = e2r(index) - 1   ! map the point from the extended mesh into the normal mesh
     dens = 0.0d0
     do imu = 1, num_orb(in1)
      dens = dens + rhoA(imu,iatom)*psi1(imu)*psi1(imu)
      renorm = renorm + dvol
      !hit = hit + 1
      !psi2 = psi2 + psi1(imu)*psi1(imu)
     enddo ! do imu
    qtot = qtot + dens*dvol
    !rhoG0(ind) = dens + rhoG0(ind)
    ewfaux(ind) = ewfaux(ind) + dens * f_mul 
    end do ! do imesh
   end do ! do iatom

   write (*,*)   "Qtot ", qtot, "renorm ", renorm, "dvol ", dvol
   return
 end subroutine project_dens0

! ===========================================================================
! ===========================================================================

 subroutine renorm_dens(ewfaux)
    use grid
    use configuration
    !use dimensions
    use interactions
    !use neighbor_map
    !use density
    use charges
    !use options
    implicit none
 ! ==== Argument Declaration and Description
    real, dimension(nrm), intent (out) :: ewfaux
    real dens,qtot,renorm
    integer i,in1,issh,iatom
 ! ==== Body
   !dens = 4.0d0*3.141592653589793238462643*Rc_max**3/3.0d0
   !write (*,'(a,5f14.7)') 'SPHERE =',dvol,renorm, dens, renorm-dens,renorm/dens
   dens = 0.0d0
   do i = 1,nrm
    dens = dens + ewfaux(i)*dvol
   enddo
   write (*,*) ' -- Total atomic density before renormalization =',dens
   qtot = 0.0d0
   do iatom = 1, natoms
     in1 = imass(iatom)
     do issh = 1, nssh(in1)
       qtot = qtot + Qneutral(issh,in1)
     end do
   end do
   renorm = qtot/dens  ! the renormalization factor
   dens = 0.0d0
   do i =1,nrm
    ewfaux(i) = ewfaux(i)*renorm
    dens      = dens + ewfaux(i)*dvol
   enddo
   write (*,*) ' -- Total atomic density after renormalization =',dens
 end subroutine renorm_dens

! ===========================================================================
! ===========================================================================

 subroutine initdenmat0( )
    use density
    use configuration
    use neighbor_map
    use interactions
    use charges
    !use options
    !use grid
    implicit none
 
 ! ======= variables
    integer imu
    !integer inu
    integer i
    integer issh
    integer iatom
    integer jatom
    integer mbeta
    integer matom
    integer in1
    !integer in2
    integer ineigh
    !integer a1
    !integer a2
    !integer a3
 
    real qmu
    logical isfile
 ! ========== Body
    rhoA = 0.0d0
 ! Loop over all atoms iatom in the unit cell
    do iatom = 1, natoms
       in1 = imass(iatom)
 ! find atom intself in list of neighbors
       matom = -99
       do ineigh = 1, neighn(iatom)
          mbeta = neigh_b(ineigh,iatom)
          jatom = neigh_j(ineigh,iatom)
          if (iatom .eq. jatom .and. mbeta .eq. 0) matom = ineigh
       end do
       imu = 1
       do issh = 1, nssh(in1)
 ! evaluate neutral charge per orbital
          qmu = Qneutral(issh,in1) / real(2*lssh(issh,in1)+1)
          do i = 1, (2*lssh(issh,in1)+1)
 ! set diagonal part of density matrix
             rhoA(imu,iatom) = qmu
             !rho(imu,imu,matom,iatom) = qmu
             imu = imu + 1
          enddo ! do i
       enddo ! do imu
    end do ! do iatom
    return
  end subroutine initdenmat0