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


! interpolate_2d.f90
! Program Description
! ===========================================================================
!       This routine is a two-dimensional interpolater on a 4x4 sub-grid.
!
! ===========================================================================
! Code rewritten by:
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
! ===========================================================================
        subroutine interpolate_2d (xin, yin, iforce, nx, ny, hx, hy, xintegral, Q_L, dQ_Ldx, dQ_Ldy)
                use dimensions
                use interactions
                use constants_fireball
                use timing
                implicit none
         
        ! Argument Declaration and Description
        ! ===========================================================================
        ! Input
                integer, intent (in) :: iforce
        
        ! Number of points in 3-center integrals
                integer, intent(in) :: nx
                integer, intent(in) :: ny
        
                real, intent (in) :: xin
                real, intent (in) :: yin
        
        ! Difference between points along x and y
                real, intent (in) :: hx
                real, intent (in) :: hy
        
                real, intent (in), dimension (numXmax, numYmax) :: xintegral
         
        ! Output
                real, intent (out) :: Q_L    ! the contibutions for a matrix element
                real, intent (out) :: dQ_Ldx ! d/dx Q_L (Computed only if iforce = 1)
                real, intent (out) :: dQ_Ldy ! d/dy Q_L (Computed only if iforce = 1)
        


        ! 1   bcna
        ! 2   xc3c
        ! 3   den3
        ! 4   den3S

        ! Local Parameters and Data Declaration
        ! ===========================================================================
        ! If you think that lower order interpolation for slowly changing parts
        ! of the surface are bad, set these two to zero.
        !        real, parameter :: tiny = 1.0d-5
        !        real, parameter :: small= 1.0d-4
        
        ! Local Variable Declaration and Description
        ! ===========================================================================
                integer imidx, imidy
                integer k, ik
        
                real px, py, inv_hx, inv_hy
                real, parameter :: xmin = 0
                real, parameter :: ymin = 0
        
                real f1m1, f0p3, f1p3, f2p1
                !real f1m2, f1m3, f0p6, f1p6
                real bb0,bb1,bb2,bb3
        ! -1 to 2 since we do third order
                !real, dimension (4,4) :: fun
                real, dimension (4)   :: g,gp

        
        ! Procedure
        ! ===========================================================================
        ! We assume xmin = ymin = 0.0d0 always in our interpolations.
        ! We assume that x<xmin and x>xmax has been checked for
        ! ===========================================================================
        ! We need to find what point of the grid to use
        
                ncall_interpolate_2d                         = ncall_interpolate_2d + 1
                ncall_interpolate_2d_inter(interaction_glob) = ncall_interpolate_2d_inter(interaction_glob) + 1 

        
                inv_hx = 1/hx
                inv_hy = 1/hy
        
                imidx = int((xin - xmin)*inv_hx)
                if (imidx .lt. 1) then
                  imidx = 1
                else if (imidx .gt. nx - 1) then
                  imidx = nx - 1
                end if
        
                imidy = int((yin - ymin)*inv_hy)
                if (imidy .lt. 1) then
                  imidy = 1
                else if (imidy .gt. ny - 1) then
                  imidy = ny - 1
                end if
        
                px = xin*inv_hx - imidx
                py = yin*inv_hy - imidy
         
                 do k = 1, 4
                  ik = k + imidx-1
                  f1m1 =   xintegral(ik,imidy     )
                  f0p3 = 3*xintegral(ik,imidy + 1 )
                  f1p3 = 3*xintegral(ik,imidy + 2 )
                  f2p1 =   xintegral(ik,imidy + 3 )
                  bb3 = -   f1m1 +   f0p3 -   f1p3 + f2p1
                  bb2 =   3*f1m1 - 2*f0p3 +   f1p3
                  bb1 = - 2*f1m1 -   f0p3 + 2*f1p3 - f2p1
                  bb0 =   2*f0p3
                  g(k)                     = ((bb3*py + bb2)*py + bb1)*py + bb0
                  if (iforce .eq. 1) gp(k) = ((3*bb3*py + 2*bb2)*py + bb1)
                 end do ! k
                 f1m1 =   g(1)
                 f0p3 = 3*g(2)
                 f1p3 = 3*g(3)
                 f2p1 =   g(4) 
                 bb3 = -   f1m1 +   f0p3 -   f1p3 + f2p1
                 bb2 =   3*f1m1 - 2*f0p3 +   f1p3
                 bb1 = - 2*f1m1 -   f0p3 + 2*f1p3 - f2p1
                 bb0 =   2*f0p3
                 Q_L = (((bb3*px + bb2)*px + bb1)*px + bb0)/36
                 if (iforce .eq. 1) then
                   dQ_Ldx = ((3*bb3*px + 2*bb2)*px + bb1)*inv_hx/36
                   f1m1 =   gp(1)
                   f0p3 = 3*gp(2)
                   f1p3 = 3*gp(3)
                   f2p1 =   gp(4) 
                   bb3 = -   f1m1 +   f0p3 -   f1p3 + f2p1
                   bb2 =   3*f1m1 - 2*f0p3 +   f1p3
                   bb1 = - 2*f1m1 -   f0p3 + 2*f1p3 - f2p1
                   bb0 =   2*f0p3
                   dQ_Ldy = (((bb3*px + bb2)*px + bb1)*px + bb0)*inv_hy/36
                end if
        
        ! Format Statements
        ! ===========================================================================
         
                return
                end