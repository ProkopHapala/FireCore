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
!       Project bands on the mesh.
!
!
!                 + X0 (iatom)
!                /|\     u1X = g1 - X0
! uX0 = X0- g0  / | \
!              /  |  + g1 (nearest grid point to iatom)
!             /   | /
!            /    |/
!           /     + Y0
!          +
!          g0 (origin of grid coords)
!
! g1 ... point where the origin of atomic mesh (sphere) is placed
!  vectors we need:
!    uX0 = X0 - g0
!    u1X = g1 - X0
!    r21 = Y0 - X0
!    u1Y = g1 - Y0 = g1 - Y0 - X0 + X0 = u1X - r21
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
 subroutine project_dens( ewfaux, f_mul )

   use configuration
   use dimensions
   use interactions
   use neighbor_map
   use grid
   use density
   use charges
   use kpoints
   implicit none

! Argument Declaration and Description
! ===========================================================================
   !integer iband0,iband1
   real, dimension (nrm), intent (out) :: ewfaux
   ! real, target, dimension (:), allocatable :: ewfaux
   real, intent (in) :: f_mul

! Local Variable Declaration and Description
! ===========================================================================

   integer iatom, jatom, ineigh
   integer imu, inu
   integer in1, in2
   integer mbeta
   integer index, index0, ind
   integer i, j, k
   integer i0, j0, k0
   integer lmu
   integer issh
   integer l
   integer imesh
   integer, dimension (3) :: ipiv
   integer, dimension (3) :: nr

   real distX, distY
   real dens

   real, dimension (3) :: r1, r2, r21
   real, dimension (3) :: u
   real, dimension (3) :: dXr, dYr
   real, dimension (3) :: X0, Y0
   real, dimension (3) :: g1
   real, dimension (3) :: u1X
   real, dimension (3) :: uX0

!   real, dimension (nspec_max):: rcutoff_max
   real, dimension (numorb_max)   :: psi1,psi2
   real, dimension (3,numorb_max) :: dpsi1,dpsi2
   real :: psiR, dpsiR
   real, dimension (  5)            :: psiL
   real, dimension (3,5)          :: dpsiL

   !real, target, dimension (:), allocatable  :: rhotmp
   real, dimension (3,3)          :: amat
   real, dimension (3,3)          :: inva
   real, dimension (:), pointer   :: pmat
   !character (len=40) filename
   !character (len=30) mssg

   !real vmin,vmax

! Procedure
! ===========================================================================

    !ewfaux = 0.0d0

    nr(1) = rm1
    nr(2) = rm2
    nr(3) = rm3
    amat = transpose(elvec)
    call inv3x3 (amat,inva)

       do iatom = 1, natoms
        in1 = imass(iatom)
        r1(:) = ratom(:,iatom)
    ! vector between the atom and initial point of the grid
        do i = 1,3
         u(i) = ratom(i,iatom) - g0(i)
        enddo ! i
    ! get n-vector
        call mult3x1 (inva,u)
    ! round coefficients to get the position of the nearest grid point g1 to the iatom X1
    ! i,j,k can be positive or negative exceeding rmX (it means not centered in the unit cell)
        i0 = nint( u(1) )
        j0 = nint( u(2) )
        k0 = nint( u(3) )
    ! find the vector u1 between the iatom X1 and the nearest point g1
        u1X(1) = u(1) - real(i0)
        u1X(2) = u(2) - real(j0)
        u1X(3) = u(3) - real(k0)
    
    ! check if the nearest grid point is located within the unit cell of the grid coords
    ! if not, let's map it within
    !i0
        if (u(1) .lt. 0.0d0) then
         i0 = i0 + rm1*(int(abs(i0/rm1)) + 1)
        else
         i0 = i0 - rm1*int(i0/rm1)
        endif
    !j0
        if (u(2) .lt. 0.0d0) then
         j0 = j0 + rm2*(int(abs(j0/rm2)) + 1)
        else
         j0 = j0 - rm2*int(j0/rm2)
        endif
    !k0
        if (u(3) .lt. 0.0d0) then
         k0 = k0 + rm3*(int(abs(k0/rm3)) + 1)
        else
         k0 = k0 - rm3*int(k0/rm3)
        endif
    ! find the coordinates of the nearest point g1 witihin the grid coords
        g1(1) = i0*elvec(1,1) + j0*elvec(2,1) + k0*elvec(3,1)
        g1(2) = i0*elvec(1,2) + j0*elvec(2,2) + k0*elvec(3,2)
        g1(3) = i0*elvec(1,3) + j0*elvec(2,3) + k0*elvec(3,3)
    ! evaluate coordinates of the iatom in the grid coords
        X0(1) = g1(1) + u1x(1)*elvec(1,1) + u1x(2)*elvec(2,1) + u1x(3)*elvec(3,1)
        X0(2) = g1(2) + u1x(1)*elvec(1,2) + u1x(2)*elvec(2,2) + u1x(3)*elvec(3,2)
        X0(3) = g1(3) + u1x(1)*elvec(1,3) + u1x(2)*elvec(2,3) + u1x(3)*elvec(3,3)
    ! vector pointing from g1 to X0
        u1X(1) = g1(1) - X0(1)
        u1X(2) = g1(2) - X0(2)
        u1X(3) = g1(3) - X0(3)
    ! save iatom coord within the grid unit cell
        ratom2g(:,iatom) = X0(:)
    ! find index of the gX point within the extended mesh
        index0 = 1 + (i0+emx1) + em1*(j0+emx2) + em1*em2*(k0+emx3)
    ! Loop over the neighbors
        do ineigh = 1, neighn(iatom)
         jatom = neigh_j(ineigh,iatom)
          mbeta = neigh_b(ineigh,iatom)
          r2(:) = ratom(:,jatom) + xl(:,mbeta)
          in2 = imass(jatom)
          do i = 1,3
           r21(i) = r2(i) - r1(i)
          enddo
    ! Loop over points in the atomic mesh gP
          do imesh = 1, nam
           index = index0 + am2rc(imesh)  ! restore index of the given mesh point gP within the extended mesh
           do i = 1,3
            dXr(i) = ram2rc(i,imesh) + u1X(i)   ! evaluate the vector between the iatom and the mesh point gP
           enddo
           do i = 1,3
            dYr(i) = dXr(i) - r21(i)   ! evaluate the vector between the jatom and the mesh point gP
           enddo
           distY = sqrt(dYr(1)**2 + dYr(2)**2 + dYr(3)**2)  ! distance between the mesh point and jatom
           if (distY .lt. Rc_max) then
            distX = sqrt(dXr(1)**2 + dXr(2)**2 + dXr(3)**2)  ! distance between the mesh point and iatom
            psi1 = 0.0d0
    ! restore wavefunctions of iatom
            imu = 1
            do issh = 1,nssh(in1)
             call getpsi(in1,issh,distX,psiR,dpsiR)   ! get radial part of wf.
             l = lssh(issh,in1)
             call getYlm(l,dXr,psiL,dpsiL)            ! get spherical harmonics part of wf.
             do lmu = 1, (2*l+1)
              psi1(imu) = psi1(imu) + psiL(lmu)*psiR
              imu = imu + 1
             enddo ! do lmu
            enddo ! do issh
            psi2 = 0.0d0
            imu = 1
            do issh = 1,nssh(in2)      ! restore wavefunctions of jatom
             call getpsi(in2,issh,distY,psiR,dpsiR)   ! get radial part of wf.
             l = lssh(issh,in2)
             call getYlm(l,dYr,psiL,dpsiL)            ! get spherical harmonics part of wf.
             do lmu = 1, (2*l+1)
              psi2(imu) = psi2(imu) + psiL(lmu)*psiR   
              imu = imu + 1
             enddo ! do lmu
            enddo ! do issh
            ind = e2r(index) - 1  ! map the point from the extended mesh into the normal mesh
            dens = 0.0d0
            do inu = 1, num_orb(in1)
             do imu = 1, num_orb(in2)
              dens = dens + rho(inu,imu,ineigh,iatom)*psi1(inu)*psi2(imu)
             enddo ! do inu
            enddo ! do imu
            ewfaux(ind) = ewfaux(ind) + dens * f_mul ! store variation of density at given point
           endif ! if (Rc_max)
          end do ! do imesh
        end do ! do ineigh
       end do ! do iatom

    !write (*,*) "project_orb vmin, vmax ", vmin, vmax
   return
 end subroutine project_dens

