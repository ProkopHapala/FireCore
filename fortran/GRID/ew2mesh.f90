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
 subroutine ew2mesh ! (icluster)
   use options, only: icluster
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
! Input
!   integer, intent (in) :: icluster

!Output


! Local Parameters and Data Declaration
! ===========================================================================
   interface
    subroutine writeout_xsf (xsfname, message, aa)
     real, dimension (:), pointer, intent (in) :: aa
     character (len=40) xsfname
     character (len=30) message
    end
  end interface

! Local Variable Declaration and Description
! ===========================================================================
   integer iatom
   integer ikpoint
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
   integer lmu
   integer issh
   integer l
   integer imesh
   integer job
   integer, dimension (3) :: ipiv
   integer, dimension (3) :: nr
   integer ii
   integer nnu
   integer mmu
   integer file
   integer iband
   integer nmax

   real distX
   real distY
   real dens
!   real densmax
   real dot
   real gutr

   real, dimension (3) :: r1
   real, dimension (3) :: r2
   real, dimension (3) :: r21
   real, dimension (3) :: u
   real, dimension (3) :: dXr
   real, dimension (3) :: dYr
   real, dimension (3) :: X0
   real, dimension (3) :: Y0
   real, dimension (3) :: g1
   real, dimension (3) :: u1X
   real, dimension (3) :: uX0



!   real, dimension (nspec_max):: rcutoff_max
   real, dimension (numorb_max)   :: psi1
   real, dimension (3,numorb_max) :: dpsi1
   real, dimension (numorb_max)   :: psi2
   real, dimension (3,numorb_max) :: dpsi2
   real :: psiR
   real :: dpsiR
   real, dimension (5)            :: psiL
   real, dimension (3,5)          :: dpsiL

!   real, dimension (3,natoms)     :: ratom2g
   real, dimension (3,3)          :: lmat
   real, dimension (3,3)          :: invl

   real,    target, dimension (:), allocatable :: ewfaux
   integer nbandsin
   integer, dimension         (:), allocatable :: pbandsin
   integer, dimension         (:), allocatable :: pkpointsin

   complex phase
   complex step1
   complex step2
   complex ai

   character(40)   :: namewf
   character(4)    :: name
   real, dimension (:), pointer   :: pmat
   character (len=30) mssg

! Procedure
! ===========================================================================

   !write (*,*) " HEY ew2mesh here, iewform = ", iewform
   !if (  iewform .eq. 5)                            call  ew2mesh_fourier (icluster)
   !if (  iewform .eq. 6)                            call  ew2mesh_kscan (icluster)
   !if (  iewform .eq. 7)                            call  ew2mesh_ARPES (icluster)
   !if (( iewform .eq. 3) .OR. ((iewform .eq. 4)))   call  ew2mesh_gamma (icluster)

! reset variables
   drhoG = 0.0d0
   job = 0

! allocate aux arrays
   allocate ( ewfaux(0:nrm-1))

   ai = cmplx(0.0d0, 1.0d0)

! set nr(:)
   nr(1) = rm1
   nr(2) = rm2
   nr(3) = rm3

   if (iewform .eq. 2) then
    ! count bands within energy window
    nbandsin = 0
    do iband = 1, norbitals_new
      do ikpoint = 1, nkpoints
         if ((ewfewin_min .le. eigen_k(iband,ikpoint)) .and. (ewfewin_max .ge. eigen_k(iband,ikpoint))) then
         	nbandsin = nbandsin + 1
         endif
      enddo  ! do ikpoint
    enddo ! do iband
    if (nbandsin .eq. 0) then
     write (*,*) '   +++++++++++++++++++++++++++++++++++++++++++++++'
     write (*,*) '            NO bands within the interval '
     write (*,*) '              skip the band projection '
     write (*,*) '   +++++++++++++++++++++++++++++++++++++++++++++++'
     return
    endif
    allocate ( pbandsin   (nbandsin+1) )
    allocate ( pkpointsin (nbandsin+1) )
    pbandsin   (nbandsin+1) = -1
! find bands within energy window ( 2nd pass )
    ii = 1
    do iband = 1, norbitals_new
      do ikpoint = 1, nkpoints
         if ((ewfewin_min .le. eigen_k(iband,ikpoint)) .and. (ewfewin_max .ge. eigen_k(iband,ikpoint))) then
      		pbandsin  (ii) = iband
      		pkpointsin(ii) = ikpoint
            ii = ii + 1
         	write (*,'(A,i6,A,i6,A,f10.5)') 'Selected band ',iband,'  kpoint ',ikpoint,' E = ',eigen_k(iband,ikpoint) 
         endif
      enddo  ! do ikpoint
    enddo ! do iband
   else ! (iewform .eq. 2)
     nbandsin = npbands*nkpoints
     allocate ( pbandsin   (nbandsin+1) )
     allocate ( pkpointsin (nbandsin+1) )
     pbandsin   (nbandsin+1) = -1
     ii = 1
     do iband = 1,npbands
     	do ikpoint = 1,nkpoints
      	pbandsin  (ii) = pbands(iband)
      	pkpointsin(ii) = ikpoint
        ii = ii + 1
		end do  
     end do
   endif ! (iewform .eq. 2)

! ===============================
! =====     Projection     ======
! ===============================
   do ii = 1, nbandsin
    iband   = pbandsin   (ii)
    ikpoint = pkpointsin (ii)

    call project_orb(iband,ewfaux)

! ===========================
! write out eigenfunctions into bandplotXXX.xsf file (format of xcrysden visual code)
!  for details see www.xcrysden.org

    if ((iewform .eq. 1) .and. ( iband .ne. pbandsin(ii+1) )) then
     file = 100 + iband
     write (name,'(i4.4)') iband
     namewf = 'bandplot_'//name//'.xsf'
     write (*,*) '  writting down band no.',iband,' into the file ',namewf
     pmat => ewfaux
     mssg = 'density_3D'
     call writeout_xsf (namewf, mssg, pmat)
     ewfaux = 0.0d0
    end if ! iewform =1

   enddo !do ii (iband)

   if (iewform .eq. 2) then
    namewf = 'bandplot_EW.xsf'
    pmat => ewfaux
    mssg = 'density_3D'
    call writeout_xsf (namewf, mssg, pmat)
   endif ! iewform = 2

   deallocate (ewfaux)
   deallocate (pbandsin)
   

! Format Statements
! ===========================================================================
100 format (3f14.6,4x,f14.6)
200 format (e14.6)

   return
 end subroutine ew2mesh

