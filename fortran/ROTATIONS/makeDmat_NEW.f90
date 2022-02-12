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
         
                real, intent(in) :: dmat (5, 5)
                real, intent(in) :: dmatm (3, numorb_max, numorb_max)
                real, intent(in) :: ddmat (3, 5, 5)
                real, intent(in) :: matm (numorb_max, numorb_max)
                real, intent(in) :: pmat (3, 3)
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
         
                real dleft (3, 5, 5)
                real dright (3, 5, 5)
                real left (5, 5)
                real right (5, 5)

                ! --- Version 1
                ! integer nl2,nl1, nm1, nm2
                ! real :: left_i,right_i,matm_i
                ! real :: dright_i(3)
                ! real :: dmatm_i(3)
                ! real :: wterm (3, 5, 5)

                ! --- Version 2
                integer nl2,nl1
                real matm_n (  5,5)
                real dmatm_n(3,5,5)
                real dli(3)
                real wti(3)
                real mi,li,ri

                ! --- Version 3
                ! integer nl2,nl1
                ! integer nm1, nm2
                ! real mi,li,ri
                ! real wti(3)
         
        ! Procedure
        ! ===========================================================================
        
                ncall_makeDmat=ncall_makeDmat+1 
        
                n1 = 0
                do issh = 1, nssh(in1)
                 l1 = lssh(issh,in1)
                 call chooserd (l1, ddmat, dpmat, dleft)
                 call chooser  (l1,  dmat,  pmat,  left)
         
                 n2 = 0
                 do jssh = 1, nssh(in2)
                  l2 = lssh(jssh,in2)
                  call chooserd (l2, ddmat, dpmat, dright)
                  call chooser  (l2,  dmat,  pmat,  right)
         
                  nl1 = 2*l1 + 1
                  nl2 = 2*l2 + 1                  
                  
        !  --- Vesrion 1  
        !         wterm(:,:nl1,:nl2) =0
        !         do m2 = 1, nl2
        !         nm2 = n2+m2     
        !         do k2 = 1, nl2
        !          right_i  = right (  k2,m2)
        !          dright_i = dright(:,k2,m2)
        !          do m1 = 1, nl1 
        !           nm1 = n1+m1  
        !           matm_i     = matm  (  nm1,nm2)
        !           dmatm_i(:) = dmatm (:,nm1,nm2)
        !           do k1 = 1, nl1
        !            left_i = left(k1,m1)
        !            wterm(:,k1,k2) = wterm(:,k1,k2)  &
        !   &           + dleft(:,k1,m1) *( right_i *matm_i   )  &
        !   &           + dright_i       *( left_i  *matm_i   )  &
        !   &           + dmatm_i        *( left_i  *right_i  )
        !           end do ! k1
        !          end do ! m1
        !         end do ! k2
        !        end do ! m2
        !        term(:, (n1+1):(n1+nl1), (n2+1):(n2+nl2) ) = wterm(:,:nl1,:nl2) 

        !  --- Vesrion 2
                 matm_n(  :nl1,:nl2) =  matm(    (n1+1):(n1+nl1), (n2+1):(n2+nl2) ) 
                dmatm_n(:,:nl1,:nl2) = dmatm( :, (n1+1):(n1+nl1), (n2+1):(n2+nl2) )  
                do k1 = 1, nl1 
                do k2 = 1, nl2
                wti(:) = 0  
                do m1 = 1, nl1
                li  =  left(  k1,m1)
                dli = dleft(:,k1,m1)          
                do m2 = 1, nl2
                   mi = matm_n(m1,m2)
                   ri = right (k2,m2)
                   wti(:) = wti(:)  &
                   &   + dli              *( ri *mi )    &
                   &   + dright (:,k2,m2) *( li *mi )    &
                   &   + dmatm_n(:,m1,m2) *( li *ri )
                  end do ! m2
                 end do ! m1
                 term(:,n1+k1,n2+k2) = wti(:)
                end do ! k2
               end do ! k1

        !   --- Vesrion 3
        !         do k1 = 1, nl1 
        !         do k2 = 1, nl2
        !         wti(:) = 0  
        !         do m1 = 1, nl1
        !         nm1 = n1+m1 
        !         li  =  left(  k1,m1)        
        !         do m2 = 1, nl2
        !            nm2 = n2+m2 
        !            mi = matm (nm1,nm2)
        !            ri = right(k2,m2)
        !            wti(:) = wti(:)  &
        !            &   + dleft (:,k1,m1)  *( ri *mi )    &
        !            &   + dright(:,k2,m2)  *( li *mi )    &
        !            &   + dmatm (:,nm1,nm2)*( li *ri )
        !           end do ! m2
        !          end do ! m1
        !          term(:,n1+k1,n2+k2) = wti(:)
        !         end do ! k2
        !        end do ! k1

              
                  n2 = n2 + nl2
                 end do
                 n1 = n1 + nl1
                end do
         
        ! Format Statements
        ! ===========================================================================
         
                return
                end
        