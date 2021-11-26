! copyright info:
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


! checksum_options.f90
! Program Description
! ===========================================================================
! A subroutine checks consistency of input settings
! ===========================================================================
! Code rewritten by:
! P. Jelinek
! Department of Thin Films
! Institute of Physics AS CR
! Cukrovarnicka 10
! Prague 6, CZ-162
! FAX +420-2-33343184
! Office telephone  +420-2-20318528
! email: jelinekp@fzu.cz
!
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine checksum_options ()

        !use outputs
        use options
        use configuration
        !use MD

        implicit none

! Argument Declaration and Description
! ===========================================================================
! Input


! Local Parameters and Data Declaration
! ===========================================================================

! Local Variable Declaration and Description
! ===========================================================================
        integer checksum


! Procedure
! ===========================================================================


! CHECKSUMS
        !checksum = idogs + iharris + ihubbard + iks
        checksum = idogs + iharris 
        if (checksum .ne. 1) then
         write (*,*) '  '
         write (*,*) ' Input Options are non-consistent !!'
         write (*,*) ' use DOGS or Ext. Hubbard or Harris or KS option'
         write (*,*) ' MUST STOP HERE !!'
         stop
        endif
! assign a value to the variable itheory
        if (iharris .eq. 1) itheory = 0
        if (idogs .eq. 1) itheory = 1
        !if (ihubbard .eq. 1) itheory = 2
!        if (iwrtdipole .eq. 1) igrid = 1
!        if (iks .eq. 1) itheory = 3
!        if (iks .eq. 1) igrid = 1
! switch of iwrtden if iks (it does not make sense to do thinks twice;
! more there is problem with renormalizaiton)
!        if (iks .eq. 1) iwrtden = 0
!        if (iwrtden .eq. 1) igrid = 1
!        if (iwrtewf .eq. 1) igrid = 1
! check if Mulliken charges are selected, if not select them
!        if (iks .eq. 1 .and. iqout .ne. 2) then
!         write (*,*) ' ----- WARNINNG ----'
!         write (*,*) ' KS theory requires Mulliken charges'
!         write (*,*) ' Mulliken charges switched on'
!         iqout = 2
!        endif
         
!        if (ixczw .eq. 1) then
!            ihorsfield = 0
!            imcweda = 0
!            igsn = 0
!            iks = 0
!            iqout = 4
!            itheory_xc = 4
!        end if

        checksum = imcweda + igsn 
        !checksum = ihorsfield + imcweda + igsn + iks+ixczw
        if (checksum .ne. 1) then
         write (*,*) '  '
         write (*,*) ' Input Options are non-consistent !!'
         write (*,*) ' use Horsfield or McWEDA or generalized SN method'
         write (*,*) ' MUST STOP HERE !!'
         stop
        endif
! assign a value to the variable itheory_xc
!        if (ihorsfield .eq. 1) itheory_xc = 0
        if (igsn .eq. 1) itheory_xc = 1
        if (imcweda .eq. 1) itheory_xc = 2

        write (*,100)

! Format Statements
! ===========================================================================
100     format (2x, 70('='))

        return
      end subroutine checksum_options
