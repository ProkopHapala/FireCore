! copyright info:
!
!                             @Copyright 2001
!                           Fireball Committee
! Brigham Young University - James P. Lewis, Chair
! Arizona State University - Otto F. Sankey
! University of Regensburg - Juergen Fritsch
! Universidad de Madrid - Jose Ortega

! Other contributors, past and present:
! Auburn University - Jian Jun Dong
! Arizona State University - Gary B. Adams
! Arizona State University - Kevin Schmidt
! Arizona State University - John Tomfohr
! Lawrence Livermore National Laboratory - Kurt Glaesemann
! Motorola, Physical Sciences Research Labs - Alex Demkov
! Motorola, Physical Sciences Research Labs - Jun Wang
! Ohio University - Dave Drabold

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


! trescentros.f90
! Program Description
! ===========================================================================
!       This subroutine calculates the (three-center) matrix elements (mu,nu).
! It calculates the matrix elements of a bondcharge with a NA interaction,
! using data which was stored in - threecint(ithet,nz,nx,kbc,jxx,ind1,ind2,ind3)
! This subroutine calculates the matrix elements molecular and crystal
! coordinates of a bondcharge with the neutral atom elsewhere (i.e. neither an
! "ontop" nor an atom-neutral atom).
!
!
!           X  ^                                * (NA)
!              !                             +
! X-Z PLANE    !                          +
!              !                      +   THETA               -----> Z
!              O ---------------------------------------------O
!            1 END OF BC                                OTHER END OF BC
!               MU                                          NU
!
! The KEY for isorp for x/c 3-center. Evaluate the spherically averaged wave
! functions for all orbitals of the two atoms.
! (x y z) = charge on (1 2 3).       1=Left  2=right  3=NA.
!
! isopr = 0    Neutral  (0 0 0)
!         1             (- 0 0)
!         2             (+ 0 0)
!         3             (0 - 0)
!         4             (0 + 0)
!         5             (0 0 -)
!         6             (0 0 +)
!
!
! ymin,ymax,numy: bond charge  distances grid
! xmin,xmax,numx: neutral atom distances grid
!
! Here are the data file characteristics:
!
! ===========================================================================
! Code originally written by Jose Ortega and Juergen Fritsch
 
! Code rewritten by:
! James P. Lewis
! Department of Physics and Astronomy
! Brigham Young University
! N233 ESC P.O. Box 24658
! Provo, UT 84602-4658
! FAX (801) 422-2265
! Office Telephone (801) 422-7444
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine trescentros_ME(interaction, isorp, maxtype, in1, in2, indna,  x, y, cost, eps, bcnax, nspecies)
        use dimensions
        use interactions
        use integrals
        use timing
        implicit none
 
! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: in1
        integer, intent (in) :: in2
        integer, intent (in) :: indna
        integer, intent (in) :: interaction
        integer, intent (in) :: isorp
        integer, intent (in) :: maxtype
        integer, intent (in) :: nspecies

        real, intent (in) :: cost
        real, intent (in) :: x
        real, intent (in) :: y
 
        real, intent (in), dimension (3, 3) :: eps

! Output
        real, intent (out), dimension (numorb_max, numorb_max) :: bcnax
 
! Local Variable Declaration and Description
! ===========================================================================
        integer imu
        integer iME
        integer index
        integer inu
        integer kforce
        integer nl

        real argument
        real cost2
        !real dQ_Ldx(5)
        !real dQ_Ldy(5)
        !real Q_L(5)
        real sint
        integer nx, ny
        real xxmax, yymax
        real hx, hy
 
        real, dimension (0:ntheta - 1, ME3c_max) :: dxbcnalist
        real, dimension (0:ntheta - 1, ME3c_max) :: dybcnalist
        real, dimension (0:ntheta - 1, ME3c_max) :: bcnalist
        real, dimension (numorb_max, numorb_max) :: bcnam
        real, dimension (ME3c_max) :: hlist
        real, dimension (0:ntheta - 1) :: p
 
! Procedure
! ===========================================================================
        ncall_trescentros=ncall_trescentros+1

        kforce = 0

! Now interpolate.
! This subroutine calls the subroutine intrp1d as needed to find the value
! of the matrix elements for any given atomic configuration.
        index = icon3c(in1,in2,indna)
        if (interaction .eq. 1) then
            call t_data3c_interpolate_2d(bcna_t(isorp,index), x,y,kforce  bcnalist, dxbcnalist, dybcnalist )
        else if (interaction .eq. 2) then
            call t_data3c_interpolate_2d(xc3c_t(isorp,index), x,y,kforce, bcnalist, dxbcnalist, dybcnalist ) 
        else if (interaction .eq. 3) then
            call t_data3c_interpolate_2d(den3_t(isorp,index), x,y,kforce, bcnalist, dxbcnalist, dybcnalist )
        end if

! Now calculate the matrix elements in molecular coordinates.
! The variable cost is passed from the assembler.
        cost2 = cost**2
        argument = 1.0d0 - cost2
        if (argument .lt. 1.0d-5) argument = 1.0d-5
        sint = sqrt(argument) 

! Set up Legendre polys in cos(theta).
        if (ntheta .ne. 5) then
          write(*,*) ' ntheta must be 5, but it is ',ntheta
          stop     
        end if
        p(0) = 1.0d0
        p(1) = cost
        p(2) = (3.0d0*cost2 - 1.0d0)/2.0d0
        p(3) = (5.0d0*cost2*cost - 3.0d0*cost)/2.0d0
        p(4) = (35.0d0*cost2*cost2 - 30.0d0*cost2 + 3.0d0)/8.0d0
 
! Multiply bcnalist pieces by the appropriate Legendre polynomials and
! form hlist.
        do iME = 1, index_max3c(in1,in2)
         hlist(iME) = 0.0d0
         do nl = 0, ntheta - 1
          hlist(iME) = hlist(iME) + p(nl)*bcnalist(nl,iME)
         end do
         if (mvalue(iME,in1,in2) .eq. 1) then
          hlist(iME) = hlist(iME)*sint
         end if
        end do
 
! Now recover bcnam which is a two-dimensional array from hlist which
! is only one-dimensional.
        call recover_3c (in1, in2, hlist, bcnam)

! Rotate bcnam into crystal-coordinates: bcnam => bcnax
        call rotate_fb (in1, in2, eps, bcnam, bcnax)
 
! Format Statements
! ===========================================================================
 
        return
      end subroutine trescentros_ME
