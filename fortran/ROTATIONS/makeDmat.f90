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


! makeDmat.f90
! Program Description
! ===========================================================================
!       This routine assembles components of the matrix derivative.
!
! ===========================================================================
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
        subroutine makeDmat (in1, in2, matm, dmatm, dmat, pmat, ddmat, dpmat, term)
                use dimensions
                use interactions
                use timing
                implicit none
         
        ! Argument Declaration and Description
        ! ===========================================================================
        ! Input
                integer, intent(in) :: in1, in2
         
                real, intent(in) :: dmat  (5, 5)
                real, intent(in) :: dmatm (3, numorb_max, numorb_max)
                real, intent(in) :: ddmat (3, 5, 5)
                real, intent(in) :: matm  (numorb_max, numorb_max)
                real, intent(in) :: pmat  (3, 3)
                real, intent(in) :: dpmat (3, 3, 3)
         
        ! Output
                real, intent(out) :: term (3, numorb_max, numorb_max)
         
        ! Local Variable Declaration and Description
        ! ===========================================================================
                integer issh, jssh
                integer ix
                integer k1, k2
                integer l1, l2
                integer m1, m2
                integer n1, n2
                integer nl1,nl2
         
                real dleft (3, 5, 5)
                real dright(3, 5, 5)
                real left     (5, 5)
                real right    (5, 5)
         
        ! Procedure
        ! ===========================================================================
        
                ncall_makeDmat=ncall_makeDmat+1 
        
                n1 = 0
                do issh = 1, nssh(in1)
                 l1 = lssh(issh,in1)
                 call chooserd (l1, ddmat, dpmat, dleft)
                 call chooser (l1, dmat, pmat, left)
         
                 n2 = 0
                 do jssh = 1, nssh(in2)
                  l2 = lssh(jssh,in2)
                  call chooserd (l2, ddmat, dpmat, dright)
                  call chooser (l2, dmat, pmat, right)

                  nl1 = 2*l1 + 1
                  nl2 = 2*l2 + 1   

         
        ! Initialize
                  do k2 = 1, nl2
                   do k1 = 1, nl1
                     term(:,n1+k1,n2+k2) = 0
                   end do
                  end do
         
        ! Calculate
                  do m2 = 1, nl2
                   do k2 = 1, nl2
                    do m1 = 1, nl1
                     do k1 = 1, nl1
        !               do ix = 1, 3
        !                 term(ix,n1+k1,n2+k2) = term(ix,n1+k1,n2+k2)  &
        !      &           + dleft(ix,k1,m1)*right(k2,m2)    *matm(n1+m1,n2+m2)    &
        !      &           + left(k1,m1)    *dright(ix,k2,m2)*matm(n1+m1,n2+m2)    &
        !      &           + left(k1,m1)    *right(k2,m2)    *dmatm(ix,n1+m1,n2+m2)
        !               end do

             term(:,n1+k1,n2+k2) = term(:,n1+k1,n2+k2)  &
             &           + dleft(:,k1,m1)      *( right(k2,m2) * matm(n1+m1,n2+m2) )   &
             &           + dright(:,k2,m2)     *( left(k1,m1)  * matm(n1+m1,n2+m2) )   &
             &           + dmatm(:,n1+m1,n2+m2)*( left(k1,m1)  * right(k2,m2)      )
                     end do
                    end do
                   end do
                  end do
         
                  n2 = n2 + nl2
                 end do
                 n1 = n1 + nl1
                end do
         
        ! Format Statements
        ! ===========================================================================
         
                return
                end
        