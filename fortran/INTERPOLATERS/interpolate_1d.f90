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


! interpolate_1d.f90
! Program Description
! ===========================================================================
! interpolate_1d uses interpolation to find the value of
! f(x) for any x, given an array of equally spaced points for f(x)
!
! For polynomial interpolation see Mathews and Walker, p.329
!
! If norder is negative, then use cubic splines instead
!
! If doing "superspline" then order is irrelevent
!
! ===========================================================================
! Code written by:
! Kurt R. Glaesemann
! Henry Eyring Center for Theoretical Chemistry
! Department of Chemistry
! University of Utah
! 315 S. 1400 E.
! Salt Lake City, UT 84112-0850
! FAX 801-581-4353
! Office telephone 801-585-1078
! ===========================================================================
!
! Program Declaration
! =========================================================================
        subroutine interpolate_1d (interaction, isub, in1, in2, non2c, ioption, xin, yout, dfdx)
        use dimensions
        use constants_fireball
        use integrals
        use timing
        implicit none

! Argument Declaration and Description
! ===========================================================================
        integer, intent(in) :: interaction
        integer, intent(in) :: isub
        integer, intent(in) :: in1
        integer, intent(in) :: in2
        integer, intent(in) :: non2c
        integer, intent(in) :: ioption   ! Derivative or not?

        real, intent(in)  :: xin         ! x
        real, intent(out) :: yout        ! F(x)
        real, intent(out) :: dfdx        ! dF/dx
 
! Local Parameters and Data Declaration
! ===========================================================================
! tol=tolerance (may be needed to avoid roundoff error in the calling program)
! if xin > xmax but within, say, .001% of xmax then ignore
        ! real, parameter :: tol=1.0d-5
        ! real, parameter :: e6t=.166666667d0
        ! real, parameter :: e24t=.04166666667d0
        ! real, parameter :: e5t=.2d0
        ! real, parameter :: e2t5=.4d0

! Local Variable Declaration and Description
! ===========================================================================
!               interaction  subtypes   index
!      overlap       1          0          1
!      vna on        2          0..9       2..11   (2  + isorp)
!      vna atm       3          0..9       12..21  (12 + isorp)
!      vxc on        4          0..4       22..26
!      vxc atm       5          0..4       27..31
!      xccorr        6          0..4       32..36
!      vnl atm / on  7          0          37
!
!      declared by the index field ind2c !!!!!!!
!
! 1   overlap
! 2   vna_ontopl
! 3   vna_ontopr
! 4   vna_atom  
! 5   vnl       
! 6   xc_ontop
! 7   xc_atom   
! 8   xc_corr   
! 9   dipole_z  
! 10  dipole_y  
! 11  dipole_x  
! 12  coulomb   
! 13  kinetic   
! 14  nuxc      
! 15  den_ontopl
! 16  den_ontopr
! 17  den_atom  
! 18  dnuxc_ol  
! 19  dnuxc_or  
! 20  denS_ontopl
! 21  denS_ontopr
! 22  denS_atom  
! 23  overlapS  
! -------------------------------------------------------------
        !integer i
        integer imid, jxx
        !integer ileft, iright
        !integer iprod, isum
        !integer j, jj, jxx
        !integer k
        !integer nn
        integer nnum
        !integer norder
 
        real h
        real xmax
        real, parameter :: xmin = 0.0d0
        real xxp

        ! real f0p1, f0p10, f0p2, f0p3, f0p30, f0p6, f0p8, f1m1, f1m12
        ! real f1m16, f1m2, f1m3, f1m4, f1p14, f1p16, f1p2, f1p24, f1p3
        ! real f1p4, f1p6, f2m1, f2m3, f2p1, f2p6, f2p7, f3p1, f3p2, ftp
        ! real p, pden
        ! real prod, prodd 
        ! real sum, sumx
        ! real xprod, xsumoverj

        real, dimension (5) :: bb
        real, dimension (nfofx) :: pdenom
        real, dimension (nfofx) :: xx
 
! Cubic spline variables
        ! integer iam
        ! real, dimension (0:nfofx) :: a, b, c, d
        ! real, dimension (0:nfofx) :: alpha
        ! real, dimension (0:nfofx) :: L
        ! real, dimension (0:nfofx) :: mu
        ! real, dimension (0:nfofx) :: Z

! Superspline variables
        real aaa, bbb, ccc, ddd
 
! Procedure
! ===========================================================================
        ncall_interpolate_1d                    = ncall_interpolate_1d+1
        ncall_interpolate_1d_inter(interaction) = ncall_interpolate_1d_inter(interaction) +1

        jxx =  ind2c(interaction,isub)
        nnum = numz2c(jxx,in1,in2)
        xmax = z2cmax(jxx,in1,in2)
! note : the points must be equally spaced and start at 0
        h = (xmax - xmin)/(nnum - 1)
        if (xin .gt. xmax) then
          xxp = xmax - xmin
        else if (xin .eq. xmin) then
        !  if(superspline) then
            yout = splineint_2c(1,non2c,1,jxx,in1,in2)
            if (ioption .eq. 1) dfdx = 0
        !  else
        !    yout = xintegral_2c(non2c,1,jxx,in1,in2)
        !    if (ioption .eq. 1) dfdx = 0
        !  end if  ! superspline
          return  
        else
          xxp = xin - xmin
        end if

! note : imid is the point to the left (code assumes this)
        imid = int(xxp/h) + 1
        if (imid .gt. nnum) imid=nnum ! If we have gone off of the end

! Superspline
! Cubic splines:  One big "super spline"

        aaa=splineint_2c(1,non2c,imid,jxx,in1,in2)
        bbb=splineint_2c(2,non2c,imid,jxx,in1,in2)
        ccc=splineint_2c(3,non2c,imid,jxx,in1,in2)
        ddd=splineint_2c(4,non2c,imid,jxx,in1,in2)

        xxp=xxp-(imid-1)*h

        if(ioption .eq. 1) dfdx =             bbb + xxp*( ccc*2 + xxp*ddd*3 )
                           yout = aaa + xxp*( bbb + xxp*( ccc   + xxp*ddd   ) )
!       d2fdx2=2.0d0*ccc+6.0d0*ddd*xxp

! Format Statements
! ===========================================================================

        return
        end
