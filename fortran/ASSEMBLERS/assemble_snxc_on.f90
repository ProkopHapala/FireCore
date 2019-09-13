! copyright info:
!
!                             @Copyright 2004
!                           Fireball Committee
! Brigham Young University - James P. Lewis, Chair
! Arizona State University - Otto F. Sankey
! Universidad de Madrid - Jose Ortega
! Lawrence Livermore National Laboratory - Kurt Glaesemann
! Universidad de Madrid - Pavel Jelinek
! Brigham Young University - Hao Wang

! Other contributors, past and present:
! Auburn University - Jian Jun Dong
! Arizona State University - Gary B. Adams
! Arizona State University - Kevin Schmidt
! Arizona State University - John Tomfohr
! Motorola, Physical Sciences Research Labs - Alex Demkov
! Motorola, Physical Sciences Research Labs - Jun Wang
! Ohio State University - Dave Drabold
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

! assemble_snxc_2c.f90
! Program Description
! ===========================================================================
!       This routine assembles the two-center exchange-correlation
! (on-site - atom case) for the average density approximation. 
!
! This subroutine could be easily incorporated in assemble_2c.f90 (JOM)
!
! The double-counting xc (uxcdcc) is also calculated here.
!
! ===========================================================================
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
        subroutine assemble_snxc_on (natoms, nprocs, my_proc, iordern,       &
     &                               itheory, uxcdcc)
        use charges
        use density
        use dimensions
        use forces
        use interactions
        use neighbor_map
        implicit none
 
! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: iordern
        integer, intent (in) :: itheory
        integer, intent (in) :: my_proc
        integer, intent (in) :: natoms
        integer, intent (in) :: nprocs

! Output
        real, intent (out) :: uxcdcc
 
! Local Parameters and Data Declaration
! ===========================================================================
 
! Local Variable Declaration and Description
! ===========================================================================
        integer iatom
        integer iatomstart
        integer imu
        integer in1, in3
        integer inu
        integer matom
        integer mbeta
        integer natomsp
 
        real xc
 
        real, dimension (numorb_max, numorb_max) :: bcxcx

! Procedure
! ===========================================================================
! Initialize
        vxc_1c = 0.0d0
        vxc = 0.0d0
        if (itheory .eq. 1) vxc_ca = 0.0d0
        uxcdcc = 0.0d0
  
! Determine which atoms are assigned to this processor.
        if (iordern .eq. 1) then
         natomsp = natoms/nprocs
         if (my_proc .lt. mod(natoms,nprocs)) then
          natomsp = natomsp + 1
          iatomstart = natomsp*my_proc + 1
         else
          iatomstart = (natomsp + 1)*mod(natoms,nprocs)                      &
     &                + natomsp*(my_proc - mod(natoms,nprocs)) + 1
         end if
        else
         iatomstart = 1
         natomsp = natoms
        end if

! Loop over the atoms in the central cell.
!!$omp parallel do private (matom, in1, in3, bcxcx, xc)
        do iatom = iatomstart, iatomstart - 1 + natomsp
         matom = neigh_self(iatom)
         in1 = imass(iatom)
         in3 = in1
     
! Note: these loops are both over in1 indeed!  This is because the two
! wavefunctions are located on the same site, but the potential is located
! at a different site.
         if (itheory .eq. 1) then
! Dogs
          call build_ca_snxc_on (in1, iatom, bcxcx, xc)
! double-counting xc correction
!!$omp atomic
          uxcdcc = uxcdcc + xc
          in3 = in1
          do inu = 1, num_orb(in3)
           do imu = 1, num_orb(in1)
            vxc_ca(imu,inu,matom,iatom) =                                    &
     &       vxc_ca(imu,inu,matom,iatom) + bcxcx(imu,inu)
           end do
          end do
         else
! Harris + ext Hubbard
          call build_snxc_on (in1, iatom, bcxcx, xc)
! double-counting xc correction
          uxcdcc = uxcdcc + xc
          in3 = in1
          do inu = 1, num_orb(in3)
           do imu = 1, num_orb(in1)
            vxc(imu,inu,matom,iatom) =                                  &
     &       vxc(imu,inu,matom,iatom) + bcxcx(imu,inu)
           end do
          end do
         end if
        end do ! End loop over iatom.

! Deallocate
! ===========================================================================

! Format Statements
! ===========================================================================
 
        return
        end subroutine assemble_snxc_on
