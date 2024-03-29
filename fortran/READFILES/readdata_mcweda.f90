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


! readdata_mcweda.f90
! Program Description
! ===========================================================================
!       This routine reads the data of McWeda approximation
!
! ===========================================================================
!
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine readdata_mcweda ()

        use configuration
        use constants_fireball
        use interactions 
        use integrals
        use options
        !use gaussG
        use charges
        
        use dimensions

        implicit none
 
! Argument Declaration and Description
! ===========================================================================
! Output


! Local Parameters and Data Declaration
! ===========================================================================
 
 
! Local Variable Declaration and Description
! ===========================================================================
        integer interaction
        integer interaction_start
 
        integer maxtype, mintype
        integer in1, in2, in3
        integer index
        integer isorp
! Procedure
! ===========================================================================
! ============================================================================
!                              read data files
! ============================================================================

! ****************************************************************************
! About ioff2c and ioff3c - the interactions are turned on or off according to 
! the ioff2c and ioff3c arrays which are read in from file diagnostics.input. 

! Two Center Interactions:
! Note the last three short range Ewald, long range Ewald, and coulomb
! (extended Hubbard) are not interactions from datafiles like the rest.
! However, we give the option here to turn off these interactions anyways.
!                ioff - 2c overlap
!                ioff - 2c vna_ontopl
!                ioff - 2c vna_ontopr
!                ioff - 2c vna_atom-atom
!                ioff - 2c non-local
!                ioff - 2c xc_ontop
!                ioff - 2c xc_atom-atom
!                ioff - 2c xc_correction
!                ioff - 2c z-dipole
!                ioff - 2c y-dipole
!                ioff - 2c x-dipole
!                ioff - 2c coulomb
!                ioff - 2c kinetic
!                ioff - 2c extended Hubbard
!                ioff - 2c den_ontopl
!                ioff - 2c den_ontopr
!                ioff - 2c den_atom
!                ioff - 2c denS_ontopl
!                ioff - 2c denS_ontopr
!                ioff - 2c denS_atom
!                ioff - 2c_overlapS
!                ioff - 2c coulomb (extended Hubbard)
!                ioff - 2c short range Ewald
!                ioff - 2c long range Ewald
!
! Three Center Interactions:
!                ioff - 3c neutral-atom
!                ioff - 3c exchange-correlation
!                ioff - 3c average density OLSXC
!                ioff - 3c average density OLSXC (spheric)
! ****************************************************************************

        if(verbosity.gt.0) write (*,*) 'subroutine readdata_mcweda() '

! one-center
        if (verbosity .ge. 3) write (*,*) '  '
        if (verbosity .ge. 3) write (*,*) ' Calling   for 1-center exchange-correlation '
! JOM-info : the options are passed now in MODULES/options.f90
!       call read_1c (nspecies, itheory, itheory_xc, ispin, ioff2c(7))
       call read_1c (nspecies, ioff2c(7))

! two-center
        interaction_start = 1
        do interaction = interaction_start, 5
         if (verbosity .ge. 3) write (*,*) '  '
         if (verbosity .ge. 3) write (*,*) ' Calling read_2c for 2-Center Interaction # ', interaction
         if (verbosity .ge. 3) write (*,100)
         call read_2c (interaction, nspecies, itheory, ioff2c(interaction), nzx)
        end do


! horsfield-like interaction (not needed for SNXC) 
        if (itheory_xc .ne. 1) then
         !  xc_ontop
         interaction = 6
         if (verbosity .ge. 3) write (*,*) '  '
         if (verbosity .ge. 3) write (*,*) ' Calling read_2c for 2-Center Interaction # ', interaction
         if (verbosity .ge. 3) write (*,100)
         call read_2c (interaction, nspecies, itheory, ioff2c(interaction), nzx)
         ! xc_atom-atom
         interaction = 7
         if (verbosity .ge. 3) write (*,*) '  '
         if (verbosity .ge. 3) write (*,*) ' Calling read_2c for 2-Center Interaction # ', interaction
         if (verbosity .ge. 3) write (*,100)
         call read_2c (interaction, nspecies, itheory, ioff2c(interaction), nzx)

         ! xc_correction
         interaction = 8
         if (verbosity .ge. 3) write (*,*) '  '
         if (verbosity .ge. 3) write (*,*) ' Calling read_2c for 2-Center Interaction # ', interaction
         if (verbosity .ge. 3) write (*,100)
         call read_2c (interaction, nspecies, itheory,ioff2c(interaction), nzx)

        end if
        
        if (itheory .eq. 1) then
         interaction = 9
         if (verbosity .ge. 3) write (*,*) '  '
         if (verbosity .ge. 3) write (*,*) ' Calling read_2c for 2-Center Interaction # ', interaction
         if (verbosity .ge. 3) write (*,100)
         call read_2c (interaction, nspecies, itheory, ioff2c(interaction), nzx)

        end if

        if (idipole .eq. 1) then

         if (itheory .eq. 1) then
          interaction = 10
          if (verbosity .ge. 3) write (*,*) '  '
          if (verbosity .ge. 3) write (*,*) ' Calling read_2c for 2-Center Interaction # ', interaction
          if (verbosity .ge. 3) write (*,100)
          call read_2c (interaction, nspecies, itheory, ioff2c(interaction), nzx)

         end if

         if (itheory .eq. 1) then
          interaction = 11
          if (verbosity .ge. 3) write (*,*) '  '
          if (verbosity .ge. 3) write (*,*) ' Calling read_2c for 2-Center Interaction # ', interaction
          if (verbosity .ge. 3) write (*,100)
          call read_2c (interaction, nspecies, itheory, ioff2c(interaction), nzx)

         end if

        end if ! idipole .eq. 1
 

!       interaction = 10, 11 are x and y dipole
        interaction = 12
        if (verbosity .ge. 3) write (*,*) '  '
        if (verbosity .ge. 3) write (*,*) ' Calling read_2c for 2-Center Interaction # ', interaction
        if (verbosity .ge. 3) write (*,100)
        call read_2c (interaction, nspecies, itheory, ioff2c(interaction), nzx)

        interaction = 13
        if (verbosity .ge. 3) write (*,*) '  '
        if (verbosity .ge. 3) write (*,*) ' Calling read_2c for 2-Center Interaction # ', interaction
        if (verbosity .ge. 3) write (*,100)
        call read_2c (interaction, nspecies, itheory, ioff2c(interaction), nzx)

! Spherical OLSXC exchange-correlation
        do interaction = 15, 23
         if (verbosity .ge. 3) write (*,*) '  '
         if (verbosity .ge. 3) write (*,*) ' Calling read2c for 2-Center Interaction # ', interaction
         if (verbosity .ge. 3) write (*,100)

         call read_2c (interaction, nspecies, itheory, ioff2c(interaction),  nzx)

        end do

! Do not need xintegral_2c if doing superspline (use splineint_2c)
        if (superspline) deallocate (xintegral_2c)

!        if (igauss .eq. 0) then   ! IF_DEF_GAUSS_END
         interaction = 1   ! bcna
         if (verbosity .ge. 3) write (*,*) '  '
         if (verbosity .ge. 3) write (*,*) ' Calling read3c for 3-Center Interaction # ', interaction
         if (verbosity .ge. 3) write (*,100)
         call read_3c (interaction, nspecies, ioff3c(interaction), itheory, nzx)
         interaction = 3   ! den3 (3c - OLSXC) - average density
         if (verbosity .ge. 3) write (*,*) '  '
         if (verbosity .ge. 3) write (*,*) ' Calling read3c for 3-Center Interaction # ', interaction
         if (verbosity .ge. 3) write (*,100)
         call read_3c (interaction, nspecies, ioff3c(interaction), itheory,  nzx)
         interaction = 4   ! den3 (3c - OLSXC) - spherical average density
         if (verbosity .ge. 3) write (*,*) '  '
         if (verbosity .ge. 3) write (*,*) ' Calling read3c for 3-Center Interaction # ', interaction
         if (verbosity .ge. 3) write (*,100)
         call read_3c (interaction, nspecies, ioff3c(interaction), itheory, nzx)
         if (verbosity .ge. 3) write (*,100)
! Set up some tables for the 2d interpolator
         call setterp_2d (nspecies, itheory_xc, itheory)
!        end if ! end if igauss  ! IF_DEF_GAUSS_END



        !  hx = hx_bcna(isorp,index)
        !  hy = hy_bcna(isorp,index)
        !  nx = numx3c_bcna(isorp,index)
        !  ny = numy3c_bcna(isorp,index)
        !  xxmax = x3cmax_bcna(isorp,index)
        !  yymax = y3cmax_bcna(isorp,index)
         
        !  write(*,*) "numx3c_bcna(:,:)", numx3c_bcna(:,:)
        !  write(*,*) "numy3c_bcna(:,:)", numy3c_bcna(:,:)
        !  write(*,*) "numx3c_xc3c(:,:)", numx3c_xc3c(:,:)
        !  write(*,*) "numy3c_xc3c(:,:)", numy3c_xc3c(:,:)
        !  write(*,*) "numx3c_den3(:,:)", numx3c_den3(:,:)
        !  write(*,*) "numy3c_den3(:,:)", numy3c_den3(:,:)

        !  write(*,*) "hx_bcna(:,:)", hx_bcna(:,:)
        !  write(*,*) "hy_bcna(:,:)", hy_bcna(:,:)
        !  write(*,*) "hx_xc3c(:,:)", hx_xc3c(:,:)
        !  write(*,*) "hy_xc3c(:,:)", hy_xc3c(:,:)
        !  write(*,*) "hx_den3(:,:)", hx_den3(:,:)
        !  write(*,*) "hy_den3(:,:)", hy_den3(:,:)

        !  write(*,*) "x3cmax_bcna(:,:)", x3cmax_bcna(:,:)
        !  write(*,*) "y3cmax_bcna(:,:)", y3cmax_bcna(:,:)
        !  write(*,*) "x3cmax_xc3c(:,:)", x3cmax_xc3c(:,:)
        !  write(*,*) "y3cmax_xc3c(:,:)", y3cmax_xc3c(:,:)
        !  write(*,*) "x3cmax_den3(:,:)", x3cmax_den3(:,:)
        !  write(*,*) "y3cmax_den3(:,:)", y3cmax_den3(:,:)


        !  do interaction = 1, 4
        !         do in1 = 1, nspecies
        !                 do in2 = 1, nspecies
        !                         do in3 = 1, nspecies
        !                                 index = icon3c(in1,in2,in3)
        !                                 call getIsorpRange(interaction,itheory,in3,maxtype,mintype)
        !                                 do isorp = mintype, maxtype
        !                                         write(*,*) "numxy3c", interaction,in1,in2,in3, isorp , " bcna(",numx3c_bcna(isorp,index),numy3c_bcna(isorp,index) !,") xc3c(", numx3c_xc3c(isorp,index),numy3c_xc3c(isorp,index),") den3(", numx3c_den3(isorp,index),numy3c_den3(isorp,index),")"
        !                                 end do
        !                         end do
        !                 end do
        !         end do
        ! end do

! ! IF_DEF_GAUSS
!         if (igauss .eq. 1) then
!          gfactor(0,0)  = sqrt(1.0d0/4.0d0/pi)
!          gfactor(1,-1) = sqrt(3.0d0/4.0d0/pi)
!          gfactor(1,0)  = sqrt(3.0d0/4.0d0/pi)
!          gfactor(1,1)  = sqrt(3.0d0/4.0d0/pi)
!          gfactor(2,-2) = sqrt(15.0d0/4.0d0/pi)
!          gfactor(2,-1) = sqrt(15.0d0/4.0d0/pi)
!          gfactor(2,0)  = sqrt(15.0d0/12.0d0/4./pi)
!          gfactor(2,1)  = sqrt(15.0d0/4.0d0/pi)
!          gfactor(2,2)  = sqrt(15.0d0/4.0d0/pi)/2.0d0
! ! MHL
! ! NA potential
!          allocate (gcoefficientsVNA (max_alphas, nspec_max))
!          allocate (alphaVNA (max_alphas, nspec_max))
!          allocate (nalphaVNA (nspec_max))
! ! Electron density
!          allocate (gcoefficientsN (max_alphas, 0:nsh_max, nspec_max))
!          allocate (alphaN (max_alphas, 0:nsh_max, nspec_max))
!          allocate (nalphaN (0:nsh_max, nspec_max))
! ! Wavefunction
!          allocate (gcoefficientsPSI (max_alphas, 1:nsh_max, nspec_max))
!          allocate (alphaPSI (max_alphas, 1:nsh_max, nspec_max))
!          allocate (nalphaPSI (1:nsh_max, nspec_max))
! ! Wavefunction for spherical approximation
!          allocate (gcoefficientsPSIS (max_alphas, 1:nsh_max, nspec_max))
!          allocate (alphaPSIS (max_alphas, 1:nsh_max, nspec_max))
!          allocate (nalphaPSIS (1:nsh_max, nspec_max))
! ! NA potential by shell
!          allocate (gcoefficientsVNA_SH (max_alphas, 1:nsh_max, nspec_max))
!          allocate (alphaVNA_SH (max_alphas, 1:nsh_max, nspec_max))
!          allocate (nalphaVNA_SH (1:nsh_max, nspec_max))
!          allocate (R_na (1:nsh_max, nspec_max))
!          call readgaussG (nspecies, nzx )
!         end if
! ! END_DEF_GAUSS

! ! IF_DEF_ORDERN
! Broadcast some variables to other processors
!        if (iordern .eq. 1) then
!         call ME_max_bcast
!         call readdata_ordern_init (nspecies, ioff2c)
!        end if
! ! END_DEF_ORDERN

        if( ivec_3c .gt. 0) call vectorize_3c() ! PROKOP_2022/02/12

! Deallocate Arrays
! ===========================================================================

 
! Format Statements
! ===========================================================================
100     format (2x, 70('='))


        return
        end subroutine readdata_mcweda
 
