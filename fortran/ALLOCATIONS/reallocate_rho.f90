! copyright info:
!
!                             @Copyright 2006
!                           Fireball Committee
! Brigham Young University - James P. Lewis, Chair
! Arizona State University - Otto F. Sankey
! Universidad de Madrid - Jose Ortega
! Academy of Sciences of the Czech Republic - Pavel Jelinek

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

 
! reallocate_rho.f90
! Program Description
! ===========================================================================
!       This routine reallocates the arrays which store the densities. 
!
! ===========================================================================
! Code written by:
! James P. Lewis
! Department of Physics and Astronomy
! Brigham Young University
! N233 ESC P.O. Box 24658 
! Provo, UT 841184602-4658
! FAX 801-422-2265
! Office telephone 801-422-7444
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine reallocate_rho (natoms, neigh_max, neighPP_max, itheory_xc, &
                                   igrid)
        use density
        use charges
        use interactions
        use scf
        implicit none
 
! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: itheory_xc
        integer, intent (in) :: natoms
        integer, intent (in) :: neigh_max
        integer, intent (in) :: neighPP_max
        integer, intent (in) :: igrid
 
! Local Parameters and Data Declaration
! ===========================================================================
 
! Local Variable Declaration and Description
! ===========================================================================
        integer ndim
 
! Allocate Arrays
! ===========================================================================
 
! Procedure
! ===========================================================================
        deallocate (cape)
        deallocate (rho)
        allocate (cape (numorb_max, numorb_max, neigh_max, natoms))
        allocate (rho (numorb_max, numorb_max, neigh_max, natoms))

        deallocate (rhoPP)
        allocate (rhoPP (numorb_max, numorb_max, neighPP_max**2, natoms))

        deallocate (cape_es)
        deallocate (rho_es)
        allocate (cape_es (numorb_max, numorb_max, neigh_max, natoms))
        allocate (rho_es (numorb_max, numorb_max, neigh_max, natoms))

        deallocate (rhoPP_es)
        allocate (rhoPP_es (numorb_max, numorb_max, neighPP_max**2, natoms))

        
        if (itheory_xc .eq. 1 .or. itheory_xc .eq. 2 .or. itheory_xc .eq. 4) then 
         deallocate (arho_on)
         deallocate (arhoi_on)
         deallocate (arho_off)
         deallocate (arhoij_off)
         deallocate (rho_off)
         deallocate (rhoij_off)
         deallocate (rho_on)
         deallocate (rhoi_on)
         allocate (arho_on (nsh_max, nsh_max, natoms))
         allocate (arhoi_on (nsh_max, nsh_max, natoms))
         allocate (arho_off (nsh_max, nsh_max, neigh_max, natoms))
         allocate (arhoij_off (nsh_max, nsh_max, neigh_max, natoms))
         allocate (rho_off (numorb_max, numorb_max, neigh_max, natoms))
         allocate (rhoij_off (numorb_max, numorb_max, neigh_max, natoms))
         allocate (rho_on (numorb_max, numorb_max, natoms))
         allocate (rhoi_on (numorb_max, numorb_max, natoms))
        end if                                                       
! jel-grid
        if (igrid .eq. 1) then
         deallocate (rho_in)
         deallocate (rho_out)
         ndim = numorb_max*numorb_max*neigh_max*natoms
         allocate (rho_in (ndim))
         allocate (rho_out (ndim))
         if (ialgmix .eq. 4) then
           deallocate (mwe)
           deallocate (drwe)
           allocate (drwe (ndim))
           allocate (mwe (ndim))
         endif
! save the old density matrix
         deallocate (rho_old) 
         allocate (rho_old (numorb_max, numorb_max, neigh_max, natoms))
        endif
! end jel-grid

! Pulay mixing
        if ((ialgmix .eq. 4) .and. (igrid .ne. 1)) then
          deallocate (mwe)
          deallocate (drwe)
          ndim = nsh_max*natoms
          allocate (mwe (ndim))
          allocate (drwe (ndim))
        endif

! Deallocate Arrays
! ===========================================================================
 
! Format Statements
! ===========================================================================
 
        return
        end subroutine reallocate_rho
