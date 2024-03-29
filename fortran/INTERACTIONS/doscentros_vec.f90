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


! doscentros.f90
! Program Description
! ===========================================================================
!      This subroutine calculates the (two-center) matrix elements (mu,nu).
! There used to be five different routines that did this for all of the
! two-center interactions - doscentros.f, dosxcatm.f, dosxcontop.f,
! dosenatm.f, and dosenontop.f.  These have now all been reduced to one
! routine in order to make Fireball more lean.
!
!      This routine also calculates the derivative with respect to the
! position of the orbital of the BRA.
!
! ===========================================================================
! Original code written by Jose Ortega.
!
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
        subroutine doscentros_vec (interaction, isub, iforce, in1, in2, in3,  distance, eps, deps, sx, spx)
        use dimensions
        use interactions
        use timing
        implicit none
 
! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: interaction
        integer, intent (in) :: isub
 
! The variable interaction is the type of two-center integral that is needed.
! The variable isub is the subtype of interation
! Here is information from read2c
!
!           interaction    subtype
!               1            isorp = 0, 8    average density 
!               1            0               overlap
!               2            0               vna_ontopl
!               2            isorp = 1, 8    vna_ontopl shell isorp
!               3            0               vna_ontopr
!               3            isorp = 1, 8    vna_ontopr shell isorp
!               4            0               vna_atom
!               4            isorp = 1, 8    vna_atom  shell isorp
!               5            0               non-local
!               6            0, 4            xc_ontop
!               7            0, 4            xc_atom
!               8            0, 4            xc_correction
!               9            0               z-dipole
!               10           0               y-dipole
!               11           0               x-dipole
!               12           0               coulomb
!               13           0               kinetic
!               14           0               extended-hubbard
!               15           isorp = 1, 8    density_ontopl
!               16           isorp = 1, 8    density_ontopr
!               17           isorp = 1, 8    density_atom  
!  spherical density (only doscentrosS)
!               18           isorp = 1, 4    sph dens ontopl
!               19           isorp = 1, 4    sph dens_ontopr
!               20           isorp = 1, 4    sph den_atom
!               21           isorp = 0       sph overlap
 
! Derivatives are computed if and only if iforce = 1.
        integer, intent (in) :: iforce

        integer, intent (in) :: in1
        integer, intent (in) :: in2
        integer, intent (in) :: in3

        real, intent (inout) :: distance
        real, intent (in), dimension (3, 3, 3) :: deps
        real, intent (in), dimension (3, 3) :: eps
 
! Output
! The variable sx is < phi(mu,r-r1) ! vna(r-ratm) ! phi(nu,r-r1)>
! which is a representation in molecular coordinates.
! The variable spx is [d/d(r1) atmna(mu,nu)] where atmna(mu,nu) is
! < phi(mu,r-r1) ! vna(r-ratm) ! phi(nu,r-r1)>.
        real, intent (out), dimension (   numorb_max, numorb_max) :: sx
        real, intent (out), dimension (3, numorb_max, numorb_max) :: spx
 
! Local Parameters and Data Declaration
! ===========================================================================
 
! Local Variable Declaration and Description
! ===========================================================================
        integer imu
        integer inu
        !integer index
        integer in_, nME
 
        real, dimension (3) :: eta

        real, dimension (ME2c_max) :: dslist ! output list of matrix elements
        real, dimension (ME2c_max) :: slist  ! output list of derivatives of matrix elements

        real, dimension (numorb_max,numorb_max) :: sm
        real, dimension (numorb_max,numorb_max) :: spm
        real, dimension (3,numorb_max,numorb_max) :: spmx

! Procedure
! ===========================================================================
        ncall_doscentros = ncall_doscentros+1
! For the atom case, in3  = in1, but for everything else in3 = in2.
! For the ontop case, in2 = in1 (left) or in3 (right).
! Initialize sm, scam and sx to zero. ! < n1 | n2 | n3 >
        sm(:,:) = 0
        sx(:,:) = 0
        if (iforce .eq. 1) then
            spm (:,:) = 0
            spmx(:,:,:) = 0
            spx (:,:,:) = 0
        end if !  (iforce .eq. 1) 
        nME = index_max2c(in1,in3)
        if( (interaction .eq. 2) .or. (interaction .eq. 15) .or. (interaction .eq. 18) ) then
                in_ = in3
        else
                in_ = in2
        end if 
        call interpolate_1d_vec (interaction, isub, in1, in_, nME, iforce, distance, slist, dslist )
        call recover_2c (in1, in3, slist,  sm)
        call recover_2c (in1, in3, dslist, spm)
        call rotate_fb  (in1, in3, eps, sm, sx)
 
! ------------- FORCES
! Only compute derivative if and only if iforce = 1.
        if (iforce .eq. 1) then
! As long as epsilon1 is called with sighat in the second "spot" as call epsilon1(R1,sighat,spe), then eps(ix,3) = eta(ix).
         eta(:) = eps(:,3)
         do inu = 1, num_orb(in3)
          do imu = 1, num_orb(in1)
! Note that if we are calculating the on-site matrix elements, then the derivatives should be exactly zero.  This is what Otto referred to as the ferbie test.  
! For example, for the on-site overlap, we get an identity matrix, thus surely we should never use a "derivative of the unit matrix". 
! Note the minus sign. d/dr1 = - eta * d/dd.
           if (distance .gt. 1.0d-3) then
            spmx(:,imu,inu) = - eta(:)*spm(imu,inu)
           end if
          end do
         end do
         call rotated (in1, in3, eps, deps, sm, spmx, spx)
        end if

! Format Statements
! ===========================================================================
 
        return
      end subroutine doscentros_vec
