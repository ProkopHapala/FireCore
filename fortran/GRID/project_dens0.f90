! Copyright info:
!
!                             @Copyright 2005
!                     Fireball Enterprise Center, BYU
! Brigham Young University - James P. Lewis, Chair
! Arizona State University - Otto F. Sankey
! Universidad de Madrid - Jose Ortega
! Institute of Physics, Czech Republic - Pavel Jelinek

! Other contributors, past and present:
! Auburn University - Jian Jun Dong
! Arizona State University - Gary B. Adams
! Arizona State University - Kevin Schmidt
! Arizona State University - John Tomfohr
! Lawrence Livermore National Laboratory - Kurt Glaesemann
! Motorola, Physical Sciences Research Labs - Alex Demkov
! Motorola, Physical Sciences Research Labs - Jun Wang
! Ohio University - Dave Drabold
! University of Regensburg - Juergen Fritsch
!
! fireball-qmd is a free (GPLv3) open project.

! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.


! den2mesh.f90
! Program Description
! ===========================================================================
!       Project density on the mesh.
!
! ===========================================================================
! Code written by:
! ===========================================================================
! P. Jelinek
! Department of Thin Films
! Institute of Physics AS CR
! Cukrovarnicka 10
! Prague 6, CZ-162
! FAX +420-2-33343184
! Office telephone  +420-2-20318530
! email: jelinekp@fzu.cz
!
! Program Declaration
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
   integer jatom
   integer mbeta
   integer ineigh
   integer index
   integer index0
   integer ind
   integer i, j, k
   integer i0, j0, k0
   integer im, jm, km
   integer lmu
   integer issh
   integer l
   integer imesh
   integer info
   integer hit
   integer, dimension (3) :: ipiv
   integer, dimension (3) :: nr

   real distX
   real qtot
   real renorm
   real dens
   real vna0
   real psi2
   real, dimension (3) :: dvna0

   real, dimension (3) :: r1
   real, dimension (3) :: r2
   real, dimension (3) :: u
   real, dimension (3) :: dXr
   real, dimension (3) :: g1
   real, dimension (3) :: X0
   real, dimension (3) :: u1X
   real, dimension (3) :: uX0


   real, dimension (numorb_max)   :: psi1
   real, dimension (3,numorb_max) :: dpsi1
   real :: psiR
   real :: dpsiR
   real, dimension (5)            :: psiL
   real, dimension (3,5)          :: dpsiL
   real, dimension (:), pointer   :: pmat
   character (len=40) filename
   character (len=30) mssg

!   real, dimension (3,natoms)     :: ratom2g
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

! reset variables
   !rhoG0 = 0.0d0
   !vnaG = 0.0d0
   psi2 = 0.0d0
   qtot = 0.0

! set nr(:)
   nr(1) = rm1
   nr(2) = rm2
   nr(3) = rm3

! make a copy of the elem grid lattice vector
! elvec
!  | a1x a1y a1z |
!  | a2x a2y a2z |
!  | a3x a3y a3z |
! position of atom is described as:
!  x = n1*a1x + n2*a2x + n3*a3x
! we need to solve this linear eq.
!
!  | a1x  a2x  a3x |   |n1|   |x|
!  | a1y  a2y  a3y | x |n2| = |y|
!  ! a1z  a2z  a3z |   |n3|   |z|
!
! copy and invert original elvec to get form written above
   lmat = transpose(elvec)
! inverse A: solving A*n=x -> n=A-1*x
   call inv3x3 (lmat,invl)

    write (*,*) 'elvec'
    write (*,*) elvec(1,:)
    write (*,*) elvec(2,:)
    write (*,*) elvec(3,:)
! integration checking
   renorm = 0.0d0
   hit = 0

! Loop over atoms
   do iatom = 1, natoms

    in1 = imass(iatom)
    r1(:) = ratom(:,iatom)

! vector between the iatom (not centered in the unit cell yet) and
! the initial grid point
    do i = 1,3
     u(i) = ratom(i,iatom) - g0(i)
    enddo ! i

! get n-vector
    call mult3x1 (invl,u)
! check it
    u1X = u
    call mult3x1 (lmat,u1X)
    do i = 1,3
      if (abs(u1X(i)-(ratom(i,iatom) - g0(i))) .gt. 0.000001d0 ) then
         write (*,*) ' check out'
         write (*,*) 'vec'
         write (*,*) u1X(:)
         write (*,*) 'diff = ',abs(u1X(i)-(ratom(i,iatom) - g0(i)))
         write (*,*) 'stop'
         stop
       endif
    enddo

! round coefficients to get the position of the nearest grid point g1 to the iatom X1
! i,j,k can be positive or negative exceeding rmX (it means not centered in the unit cell)
    i0 = nint( u(1) )
    j0 = nint( u(2) )
    k0 = nint( u(3) )

! find the vector u1 between the iatom X1 and the nearest point g1
    u1X(1) = u(1) - real(i0)
    u1X(2) = u(2) - real(j0)
    u1X(3) = u(3) - real(k0)
    if (u(1) .lt. 0.0d0) then
     i0 = i0 + rm1*(int(abs(i0/rm1)) + 1)
    else
     i0 = i0 - rm1*int(i0/rm1)
    endif
    if (u(2) .lt. 0.0d0) then
     j0 = j0 + rm2*(int(abs(j0/rm2)) + 1)
    else
     j0 = j0 - rm2*int(j0/rm2)
    endif
    if (u(3) .lt. 0.0d0) then
     k0 = k0 + rm3*(int(abs(k0/rm3)) + 1)
    else
     k0 = k0 - rm3*int(k0/rm3)
    endif
! find the coordinates of the nearest point g1 witihin the grid coords
    !g1(1) = i0*elvec(1,1) + j0*elvec(2,1) + k0*elvec(3,1)
    !g1(2) = i0*elvec(1,2) + j0*elvec(2,2) + k0*elvec(3,2)
    !g1(3) = i0*elvec(1,3) + j0*elvec(2,3) + k0*elvec(3,3)
    g1(:) = i0*elvec(1,:) + j0*elvec(2,:) + k0*elvec(3,:)
! evaluate coordinates of the iatom in the grid coords
    !X0(1) = g1(1) + u1x(1)*elvec(1,1) + u1x(2)*elvec(2,1) + u1x(3)*elvec(3,1)
    !X0(2) = g1(2) + u1x(1)*elvec(1,2) + u1x(2)*elvec(2,2) + u1x(3)*elvec(3,2)
    !X0(3) = g1(3) + u1x(1)*elvec(1,3) + u1x(2)*elvec(2,3) + u1x(3)*elvec(3,3)
    X0(:) = g1(:) + u1x(1)*elvec(1,:) + u1x(2)*elvec(2,:) + u1x(3)*elvec(3,:)
! vector pointing from g1 to X0
    !u1X(1) = g1(1) - X0(1)
    !u1X(2) = g1(2) - X0(2)
    !u1X(3) = g1(3) - X0(3)
    u1X(:) = g1(:) - X0(:)

    ratom2g(:,iatom) = X0(:)                          ! save iatom coord within the grid unit cell
    index0 = 1 + (i0+emx1) + em1*(j0+emx2) + em1*em2*(k0+emx3)  ! find index of the gX point within the extende mesh
! skip certain atoms to avoid the overcounting; it means all jatoms having the identification number less then iatom because those pairs have been counted previously
    ! Loop over points in the atomic mesh gP
    do imesh = 1, nam
     index = index0 + am2rc(imesh)
      dXr(:) = ram2rc(:,imesh) + u1X(:)
     distX = sqrt(dXr(1)**2 + dXr(2)**2 + dXr(3)**2)  ! distance between the mesh point and iatom
     psi1 = 0.0d0
     imu = 1
     do issh = 1,nssh(in1)
      call getpsi(in1,issh,distX,psiR,dpsiR)
      l = lssh(issh,in1)    ! angular momentum
      call getYlm(l,dXr,psiL,dpsiL)
      do lmu = 1, (2*l+1)
       psi1(imu) = psi1(imu) + psiL(lmu)*psiR  ! get spherical harmonics part of wf.
       imu = imu + 1
      enddo ! do lmu
     enddo ! do issh

     ind = e2r(index) - 1   ! map the point from the extended mesh into the normal mesh
! assemble density
     dens = 0.0d0
     do imu = 1, num_orb(in1)
      dens = dens + rhoA(imu,iatom)*psi1(imu)*psi1(imu)
      renorm = renorm + dvol
      hit = hit + 1
      psi2 = psi2 + psi1(imu)*psi1(imu)
     enddo ! do imu

    qtot = qtot + dens*dvol
    !rhoG0(ind) = dens + rhoG0(ind)
    ewfaux(ind) = ewfaux(ind) + dens * f_mul 

! ==========get vna potential
!     call getvna (in1, distX, vna0, dvna0)
!     vnaG(ind) = vnaG(ind) + vna0
! store variation of density at given point
    end do ! do imesh
   end do ! do iatom

   write (*,*)   "Qtot ", qtot, "renorm ", renorm, "dvol ", dvol
   return
 end subroutine project_dens0

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