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


! allocate_f.f90
! Program Description
! ===========================================================================
!       This routine allocates the arrays which store the derivatives of the
! interactions of the Hamiltonian matrix.
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
        !subroutine allocate_f (natoms, neigh_max, neighPP_max, numorb_max, nsh_max, itheory, itheory_xc, igauss, ivdw,  iharmonic, ibias)
        subroutine allocate_f
        use forces
        use options !, only : idipole
        use configuration
        use interactions
        use neighbor_map
        implicit none
 
! Argument Declaration and Description
! ===========================================================================
! Input
        ! integer, intent (in) :: igauss
        ! integer, intent (in) :: iharmonic
        ! integer, intent (in) :: itheory
        ! integer, intent (in) :: itheory_xc
        ! integer, intent (in) :: ivdw
        ! integer, intent (in) :: ibias
        ! integer, intent (in) :: natoms
        ! integer, intent (in) :: neigh_max
        ! integer, intent (in) :: neighPP_max
        ! integer, intent (in) :: numorb_max
        ! integer, intent (in) :: nsh_max
 
! Local Parameters and Data Declaration
! ===========================================================================
 
! Local Variable Declaration and Description
! ===========================================================================
 
! Allocate Arrays
! ===========================================================================
! Allocate derivatives of interactions
        if (allocated(sp_mat)) deallocate(sp_mat)
        if (allocated(sp_mat)) deallocate(sp_mat)
        if (allocated(sp_mat)) deallocate(sp_mat)
        if (allocated(sp_mat)) deallocate(sp_mat)
        allocate (sp_mat (3, numorb_max, numorb_max, neigh_max, natoms))
        if (allocated(spVNL)) deallocate(spVNL)
        if (allocated(spVNL)) deallocate(spVNL)
        if (allocated(spVNL)) deallocate(spVNL)
        if (allocated(spVNL)) deallocate(spVNL)
        allocate (spVNL (3, numorb_max, numorb_max, neighPP_max, natoms))
        if (allocated(tp_mat)) deallocate(tp_mat)
        if (allocated(tp_mat)) deallocate(tp_mat)
        if (allocated(tp_mat)) deallocate(tp_mat)
        if (allocated(tp_mat)) deallocate(tp_mat)
        allocate (tp_mat (3, numorb_max, numorb_max, neigh_max, natoms))

! Allocate components of the forces
        if (allocated(dusr)) deallocate(dusr)
        if (allocated(dusr)) deallocate(dusr)
        if (allocated(dusr)) deallocate(dusr)
        if (allocated(dusr)) deallocate(dusr)
        allocate (dusr (3, natoms))
        if (allocated(dxcv)) deallocate(dxcv)
        if (allocated(dxcv)) deallocate(dxcv)
        if (allocated(dxcv)) deallocate(dxcv)
        if (allocated(dxcv)) deallocate(dxcv)
        allocate (dxcv (3, natoms))
        if (allocated(fro)) deallocate(fro)
        if (allocated(fro)) deallocate(fro)
        if (allocated(fro)) deallocate(fro)
        if (allocated(fro)) deallocate(fro)
        allocate (fro (3, natoms))
        if (allocated(ft)) deallocate(ft)
        if (allocated(ft)) deallocate(ft)
        if (allocated(ft)) deallocate(ft)
        if (allocated(ft)) deallocate(ft)
        allocate (ft (3, natoms))
        if (allocated(ftot)) deallocate(ftot)
        if (allocated(ftot)) deallocate(ftot)
        if (allocated(ftot)) deallocate(ftot)
        if (allocated(ftot)) deallocate(ftot)
        allocate (ftot (3, natoms))
        if (allocated(ftotold)) deallocate(ftotold)
        if (allocated(ftotold)) deallocate(ftotold)
        if (allocated(ftotold)) deallocate(ftotold)
        if (allocated(ftotold)) deallocate(ftotold)
        allocate (ftotold (3, natoms))
        if (allocated(ftotnew)) deallocate(ftotnew)
        if (allocated(ftotnew)) deallocate(ftotnew)
        if (allocated(ftotnew)) deallocate(ftotnew)
        if (allocated(ftotnew)) deallocate(ftotnew)
        allocate (ftotnew (3, natoms))
        if (allocated(fana)) deallocate(fana)
        if (allocated(fana)) deallocate(fana)
        if (allocated(fana)) deallocate(fana)
        if (allocated(fana)) deallocate(fana)
        allocate (fana (3, neigh_max, natoms))
        if (allocated(faxc)) deallocate(faxc)
        if (allocated(faxc)) deallocate(faxc)
        if (allocated(faxc)) deallocate(faxc)
        if (allocated(faxc)) deallocate(faxc)
        allocate (faxc (3, neigh_max, natoms))
        if (allocated(f3naa)) deallocate(f3naa)
        if (allocated(f3naa)) deallocate(f3naa)
        if (allocated(f3naa)) deallocate(f3naa)
        if (allocated(f3naa)) deallocate(f3naa)
        allocate (f3naa (3, natoms))
        if (allocated(f3nab)) deallocate(f3nab)
        if (allocated(f3nab)) deallocate(f3nab)
        if (allocated(f3nab)) deallocate(f3nab)
        if (allocated(f3nab)) deallocate(f3nab)
        allocate (f3nab (3, natoms))
        if (allocated(f3nac)) deallocate(f3nac)
        if (allocated(f3nac)) deallocate(f3nac)
        if (allocated(f3nac)) deallocate(f3nac)
        if (allocated(f3nac)) deallocate(f3nac)
        allocate (f3nac (3, natoms))
        if (allocated(f3nla)) deallocate(f3nla)
        if (allocated(f3nla)) deallocate(f3nla)
        if (allocated(f3nla)) deallocate(f3nla)
        if (allocated(f3nla)) deallocate(f3nla)
        allocate (f3nla (3, natoms))
        if (allocated(f3nlb)) deallocate(f3nlb)
        if (allocated(f3nlb)) deallocate(f3nlb)
        if (allocated(f3nlb)) deallocate(f3nlb)
        if (allocated(f3nlb)) deallocate(f3nlb)
        allocate (f3nlb (3, natoms))
        if (allocated(f3nlc)) deallocate(f3nlc)
        if (allocated(f3nlc)) deallocate(f3nlc)
        if (allocated(f3nlc)) deallocate(f3nlc)
        if (allocated(f3nlc)) deallocate(f3nlc)
        allocate (f3nlc (3, natoms))
        if (allocated(f3xca)) deallocate(f3xca)
        if (allocated(f3xca)) deallocate(f3xca)
        if (allocated(f3xca)) deallocate(f3xca)
        if (allocated(f3xca)) deallocate(f3xca)
        allocate (f3xca (3, natoms))
        if (allocated(f3xcb)) deallocate(f3xcb)
        if (allocated(f3xcb)) deallocate(f3xcb)
        if (allocated(f3xcb)) deallocate(f3xcb)
        if (allocated(f3xcb)) deallocate(f3xcb)
        allocate (f3xcb (3, natoms))
        if (allocated(f3xcc)) deallocate(f3xcc)
        if (allocated(f3xcc)) deallocate(f3xcc)
        if (allocated(f3xcc)) deallocate(f3xcc)
        if (allocated(f3xcc)) deallocate(f3xcc)
        allocate (f3xcc (3, natoms))
        if (allocated(fotxc)) deallocate(fotxc)
        if (allocated(fotxc)) deallocate(fotxc)
        if (allocated(fotxc)) deallocate(fotxc)
        if (allocated(fotxc)) deallocate(fotxc)
        allocate (fotxc (3, neigh_max, natoms))
        if (allocated(fotna)) deallocate(fotna)
        if (allocated(fotna)) deallocate(fotna)
        if (allocated(fotna)) deallocate(fotna)
        if (allocated(fotna)) deallocate(fotna)
        allocate (fotna (3, neigh_max, natoms))
!        allocate (ftot_dftd3 (3, natoms))
 
! Procedure
! ===========================================================================
! NOTE: here we can allocate matrix of size only neighPP_max, because 
! ontol and atomic cases are just  2 center interaction, but 
! pure 3 center interactions is different story 
        if (allocated(fanl)) deallocate(fanl)
        if (allocated(fanl)) deallocate(fanl)
        if (allocated(fanl)) deallocate(fanl)
        if (allocated(fanl)) deallocate(fanl)
        allocate (fanl (3, neighPP_max, natoms))
        if (allocated(fotnl)) deallocate(fotnl)
        if (allocated(fotnl)) deallocate(fotnl)
        if (allocated(fotnl)) deallocate(fotnl)
        if (allocated(fotnl)) deallocate(fotnl)
        allocate (fotnl (3, neighPP_max, natoms))

! ! IF_DEF_GAUSS
!        if (igauss .eq. 1) allocate (fxcro (3, neigh_max, natoms))
! ! END_DEF_GAUSS

! Allocate components of the forces - needed for either DOGS or extended-Hubbard
        if (itheory .ne. 0) then
         if (allocated(dewald)) deallocate(dewald)
         if (allocated(dewald)) deallocate(dewald)
         if (allocated(dewald)) deallocate(dewald)
         if (allocated(dewald)) deallocate(dewald)
         allocate (dewald (3, natoms, natoms))
         if (allocated(fewald)) deallocate(fewald)
         if (allocated(fewald)) deallocate(fewald)
         if (allocated(fewald)) deallocate(fewald)
         if (allocated(fewald)) deallocate(fewald)
         allocate (fewald (3, natoms))
         if (allocated(flrew)) deallocate(flrew)
         if (allocated(flrew)) deallocate(flrew)
         if (allocated(flrew)) deallocate(flrew)
         if (allocated(flrew)) deallocate(flrew)
         allocate (flrew (3, natoms))
 !        allocate (flrew_qmmm (3, natoms))
        end if

! Allocate components of the forces - needed for DOGS   
        if (itheory .eq. 1 .or. idipole .eq. 1) then
         if (allocated(dipp)) deallocate(dipp)
         if (allocated(dipp)) deallocate(dipp)
         if (allocated(dipp)) deallocate(dipp)
         if (allocated(dipp)) deallocate(dipp)
         allocate (dipp (3, numorb_max, numorb_max, neigh_max, natoms))
! JIMM
         if (allocated(dippcm)) deallocate(dippcm)
         if (allocated(dippcm)) deallocate(dippcm)
         if (allocated(dippcm)) deallocate(dippcm)
         if (allocated(dippcm)) deallocate(dippcm)
         allocate (dippcm (3, 3, numorb_max, numorb_max))
         if (allocated(dippc)) deallocate(dippc)
         if (allocated(dippc)) deallocate(dippc)
         if (allocated(dippc)) deallocate(dippc)
         if (allocated(dippc)) deallocate(dippc)
         allocate (dippc (3, 3, numorb_max, numorb_max, neigh_max, natoms))

         if (allocated(faca)) deallocate(faca)
         if (allocated(faca)) deallocate(faca)
         if (allocated(faca)) deallocate(faca)
         if (allocated(faca)) deallocate(faca)
         allocate (faca (3, neigh_max, natoms))
         if (allocated(faxc_ca)) deallocate(faxc_ca)
         if (allocated(faxc_ca)) deallocate(faxc_ca)
         if (allocated(faxc_ca)) deallocate(faxc_ca)
         if (allocated(faxc_ca)) deallocate(faxc_ca)
         allocate (faxc_ca (3, neigh_max, natoms))
         if (allocated(f3caa)) deallocate(f3caa)
         if (allocated(f3caa)) deallocate(f3caa)
         if (allocated(f3caa)) deallocate(f3caa)
         if (allocated(f3caa)) deallocate(f3caa)
         allocate (f3caa (3, natoms))
         if (allocated(f3cab)) deallocate(f3cab)
         if (allocated(f3cab)) deallocate(f3cab)
         if (allocated(f3cab)) deallocate(f3cab)
         if (allocated(f3cab)) deallocate(f3cab)
         allocate (f3cab (3, natoms))
         if (allocated(f3cac)) deallocate(f3cac)
         if (allocated(f3cac)) deallocate(f3cac)
         if (allocated(f3cac)) deallocate(f3cac)
         if (allocated(f3cac)) deallocate(f3cac)
         allocate (f3cac (3, natoms))
         if (allocated(f3xca_ca)) deallocate(f3xca_ca)
         if (allocated(f3xca_ca)) deallocate(f3xca_ca)
         if (allocated(f3xca_ca)) deallocate(f3xca_ca)
         if (allocated(f3xca_ca)) deallocate(f3xca_ca)
         allocate (f3xca_ca (3, natoms))
         if (allocated(f3xcb_ca)) deallocate(f3xcb_ca)
         if (allocated(f3xcb_ca)) deallocate(f3xcb_ca)
         if (allocated(f3xcb_ca)) deallocate(f3xcb_ca)
         if (allocated(f3xcb_ca)) deallocate(f3xcb_ca)
         allocate (f3xcb_ca (3, natoms))
         if (allocated(f3xcc_ca)) deallocate(f3xcc_ca)
         if (allocated(f3xcc_ca)) deallocate(f3xcc_ca)
         if (allocated(f3xcc_ca)) deallocate(f3xcc_ca)
         if (allocated(f3xcc_ca)) deallocate(f3xcc_ca)
         allocate (f3xcc_ca (3, natoms))
         if (allocated(fotca)) deallocate(fotca)
         if (allocated(fotca)) deallocate(fotca)
         if (allocated(fotca)) deallocate(fotca)
         if (allocated(fotca)) deallocate(fotca)
         allocate (fotca (3, neigh_max, natoms))
         if (allocated(fotxc_ca)) deallocate(fotxc_ca)
         if (allocated(fotxc_ca)) deallocate(fotxc_ca)
         if (allocated(fotxc_ca)) deallocate(fotxc_ca)
         if (allocated(fotxc_ca)) deallocate(fotxc_ca)
         allocate (fotxc_ca (3, neigh_max, natoms))
        end if

! IF_DEF_HUBBARD
! Allocate components of the forces - needed for extended Hubbard
!        if (itheory .eq. 2) then
!         allocate (fcoulomb (3, natoms))
!         allocate (fxcnu (3, natoms))
!        end if
! END_DEF_HUBBARD

! Allocate snxc forces
        if (itheory_xc .eq. 1 .or. itheory_xc .eq. 2 .or. itheory_xc .eq. 4) then
         if (allocated(spm_mat)) deallocate(spm_mat)
         if (allocated(spm_mat)) deallocate(spm_mat)
         if (allocated(spm_mat)) deallocate(spm_mat)
         if (allocated(spm_mat)) deallocate(spm_mat)
         allocate (spm_mat     (3, nsh_max,    nsh_max,    neigh_max, natoms))
         if (allocated(arhop_off)) deallocate(arhop_off)
         if (allocated(arhop_off)) deallocate(arhop_off)
         if (allocated(arhop_off)) deallocate(arhop_off)
         if (allocated(arhop_off)) deallocate(arhop_off)
         allocate (arhop_off   (3, nsh_max,    nsh_max,    neigh_max, natoms)) 
         if (allocated(arhopij_off)) deallocate(arhopij_off)
         if (allocated(arhopij_off)) deallocate(arhopij_off)
         if (allocated(arhopij_off)) deallocate(arhopij_off)
         if (allocated(arhopij_off)) deallocate(arhopij_off)
         allocate (arhopij_off (3, nsh_max,    nsh_max,    neigh_max, natoms)) 
         if (allocated(arhop_on)) deallocate(arhop_on)
         if (allocated(arhop_on)) deallocate(arhop_on)
         if (allocated(arhop_on)) deallocate(arhop_on)
         if (allocated(arhop_on)) deallocate(arhop_on)
         allocate (arhop_on    (3, nsh_max,    nsh_max,    neigh_max, natoms)) 
         if (allocated(rhop_off)) deallocate(rhop_off)
         if (allocated(rhop_off)) deallocate(rhop_off)
         if (allocated(rhop_off)) deallocate(rhop_off)
         if (allocated(rhop_off)) deallocate(rhop_off)
         allocate (rhop_off    (3, numorb_max, numorb_max, neigh_max, natoms)) 
         if (allocated(rhopij_off)) deallocate(rhopij_off)
         if (allocated(rhopij_off)) deallocate(rhopij_off)
         if (allocated(rhopij_off)) deallocate(rhopij_off)
         if (allocated(rhopij_off)) deallocate(rhopij_off)
         allocate (rhopij_off  (3, numorb_max, numorb_max, neigh_max, natoms)) 
         if (allocated(rhop_on)) deallocate(rhop_on)
         if (allocated(rhop_on)) deallocate(rhop_on)
         if (allocated(rhop_on)) deallocate(rhop_on)
         if (allocated(rhop_on)) deallocate(rhop_on)
         allocate (rhop_on     (3, numorb_max, numorb_max, neigh_max, natoms)) 
         ! OLSXC double count corr forces
         if (allocated(dxcdcc)) deallocate(dxcdcc)
         if (allocated(dxcdcc)) deallocate(dxcdcc)
         if (allocated(dxcdcc)) deallocate(dxcdcc)
         if (allocated(dxcdcc)) deallocate(dxcdcc)
         allocate (dxcdcc (3, neigh_max, natoms))
        end if

! Allocate xczw forces (double countig correction)
!         if (itheory_xc .eq. 4) allocate (dxcdcc_zw (3, neigh_max, natoms))

! Allocate vdw forces if requested.
!        if (ivdw .eq. 1) allocate (fvdw (3, natoms))

! Allocate external field forces for thermodynamic integration if requested.
!        if (iharmonic .eq. 1) allocate (fharmonic (3, natoms)) 
 
! Allocate bias forces for bias voltage field if requested.
!        if (ibias .eq. 1) allocate (fbias (3, natoms)) 

! Deallocate Arrays
! ===========================================================================
 
! Format Statements
! ===========================================================================

        return
        end subroutine allocate_f
