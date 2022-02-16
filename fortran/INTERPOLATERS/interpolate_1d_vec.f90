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
        subroutine interpolate_1d_vec (interaction, isub, in1, in2, nME, ioption, xin, yout, dfdx)
        use dimensions
        use constants_fireball
        use integrals
        use interactions, only : ME2c_max
        use timing
        implicit none

! Argument Declaration and Description
! ===========================================================================
        integer, intent(in) :: interaction
        integer, intent(in) :: isub
        integer, intent(in) :: in1
        integer, intent(in) :: in2
        integer, intent(in) :: nME
        integer, intent(in) :: ioption   ! Derivative or not?

        real, intent(in)  :: xin            ! x
        real, intent(out) :: yout(ME2c_max) ! F(x)
        real, intent(out) :: dfdx(ME2c_max) ! dF/dx
 
! Local Parameters and Data Declaration
! ===========================================================================
        integer imid, jxx
        integer nnum
        real h, xmax
        real xxp
        real, dimension (5) :: bb
        real, dimension (nfofx) :: pdenom
        real, dimension (nfofx) :: xx
 
        real, dimension (nME) :: aaa, bbb, ccc, ddd
 
! Procedure
! ===========================================================================
        ncall_interpolate_1d_vec=ncall_interpolate_1d_vec+1

        jxx  =  ind2c(interaction,isub)
        nnum = numz2c(jxx,in1,in2)
        xmax = z2cmax(jxx,in1,in2)
! note : the points must be equally spaced and start at 0
        h = xmax/(nnum - 1)
        if (xin .gt. xmax) then
          xxp = xmax
        else if (xin .le. 0) then
          yout(:) = splineint_2c(1,:nME,1,jxx,in1,in2)
          if (ioption .eq. 1) dfdx = 0
          return  
        else
          xxp = xin
        end if

        imid = int(xxp/h) + 1         ! imid is the point to the left (code assumes this)
        if (imid .gt. nnum) imid=nnum ! If we have gone off of the end
        aaa(:)=splineint_2c(1,:nME,imid,jxx,in1,in2)
        bbb(:)=splineint_2c(2,:nME,imid,jxx,in1,in2)
        ccc(:)=splineint_2c(3,:nME,imid,jxx,in1,in2)
        ddd(:)=splineint_2c(4,:nME,imid,jxx,in1,in2)

        xxp=xxp-(imid-1)*h

        if(ioption .eq. 1) dfdx(:nME) =             bbb + xxp*( ccc*2 + xxp*ddd*3 )
                           yout(:nME) = aaa + xxp*( bbb + xxp*( ccc   + xxp*ddd   ) )
!       d2fdx2=2.0d0*ccc+6.0d0*ddd*xxp

! Format Statements
! ===========================================================================

        return
        end subroutine interpolate_1d_vec
