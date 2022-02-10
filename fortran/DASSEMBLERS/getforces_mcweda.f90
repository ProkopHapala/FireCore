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


! forces_mcweda.f90
! Program Description
! ===========================================================================
!       This routine assemble forces for McWeda
!
! ===========================================================================

! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine getforces_mcweda ()

        use options
        ! use outputs
        ! use mpi_main
        use configuration
        use forces
        !use qmmm_module, only : qmmm_struct
        use timing
        implicit none

! Argument Declaration and Description
! ===========================================================================
! Input


! Local Parameters and Data Declaration
! ===========================================================================

! Local Variable Declaration and Description
! ===========================================================================


! Procedure
! ===========================================================================

! ============================================================================
!                              Dassemble_2c
! ============================================================================
          !write (*,*) '  '
          !write (*,*) '  '
          !write (*,*) ' ***************************************************** '
          !write (*,*) ' Dassemble two-center force contributions. '
          call timer_start_i(51)
          call Dassemble_2c()
          call timer_stop_i(51)
          !write (*,*) ' Dassemble two-center PP force contributions. '
          call timer_start_i(52)
          call Dassemble_2c_PP()
          call timer_stop_i(52)

!           if (ibias .eq. 1) call Dassemble_bias (nprocs, iordern)     ! IF_DEF_BIAS_END



! Call the exchange-correlation interactions based on method chosen
! (i.e. itheory_xc).
          if (itheory_xc .eq. 1) then
           !write (*,*) ' Dassemble on-site SNXC force contributions. '
           if (itheory .eq. 1) then
             call Dassemble_ca_snxc_on()
             call Dassemble_ca_snxc_2c()
           else
             call Dassemble_snxc_on()
             call Dassemble_snxc_2c()
           endif
          end if

          if (itheory_xc .eq. 2) then
           !write (*,*) ' Dassemble on-site OSLXC force contributions. '
            
           if (itheory .eq. 1) then
             call timer_start_i(53)
             call Dassemble_ca_olsxc_on()
             call timer_stop_i(53)
             call timer_start_i(54)
             call Dassemble_ca_olsxc_2c()
             call timer_stop_i(54)
           else
             call Dassemble_olsxc_on()
             call Dassemble_olsxc_2c()
           endif
          end if

          call timer_start_i(55)
          if (itheory .eq. 1) then
           !write (*,*) ' Dassemble two-center DOGS force contributions. ' 
           !if (idipole .eq. 1) then      ! IF_DIPOLE_END
           ! call Dassemble_ca_2c_dip ()  ! IF_DIPOLE_END
           !else                          ! IF_DIPOLE_END
            call Dassemble_ca_2c     ()  
           !end if                        ! IF_DIPOLE_END
          endif
          call timer_stop_i(55)

! ===========================================================================
!                               Dassemble_3c
! ===========================================================================
          !write (*,*) '  '
          !write (*,*) ' Dassemble three-center force contributions. '
          call timer_start_i(56)
          call Dassemble_3c ()
          call timer_stop_i(56)

          !write (*,*) ' Dassemble three-center PP force contributions. '
          call timer_start_i(57)
          call Dassemble_3c_PP ()
          call timer_stop_i(57)

          if (itheory .eq. 1) then
           !write (*,*) ' Dassemble three-center DOGS force contributions. '

          call timer_start_i(58)
          !if (idipole .eq. 1) then     ! IF_DIPOLE_END
          !  call Dassemble_ca_3c_dip () ! IF_DIPOLE_END
          !else                         ! IF_DIPOLE_END
            call Dassemble_ca_3c ()
          !end if                       ! IF_DIPOLE_END
          !write (*,*) ' Dassemble three-center long-range contributions. '
          call timer_stop_i(58)  

          call timer_start_i(59)
          !if (idipole .eq. 1)          ! IF_DIPOLE_END
          !  call Dassemble_lr_dip ()    ! IF_DIPOLE_END
          !else                          ! IF_DIPOLE_END
            call Dassemble_lr     ()
          !end if                        ! IF_DIPOLE_END
          call timer_stop_i(59)
! IF_DEF_QMMM
!           if (iqmmm .eq. 1) then
!             !write (*,*) ' Dassemble three-center qm/mm contributions. '
!             if (idipole .eq. 0) call Dassemble_qmmm (nprocs, iordern)
!             if (idipole .eq. 1) call D_dip (nprocs, iordern)
!           else
!             flrew_qmmm = 0.0d0
!             qmmm_struct%dxyzcl = 0.0d0
!           end if
! END_DEF_QMMM
          end if

! Call the exchange-correlation interactions based on method chosen
          if (itheory_xc .eq. 1) then
           !write (*,*) ' Dassemble off-site SN exchange-correlation forces. '
           if (itheory .eq. 1) then
! Include gaussians
            call Dassemble_ca_snxc_3c ()
           else
            call Dassemble_snxc_3c ()
           endif

          else if(itheory_xc .eq. 2) then
           !write (*,*) ' Dassemble off-site OLS exchange-correlation forces. '
            call timer_start_i(60)
           if (itheory .eq. 1) then
            call Dassemble_ca_olsxc_3c ()
           else
            call Dassemble_olsxc_3c ()
           endif
           call timer_stop_i(60)
          end if

          !write (*,*) ' ***************************************************** '

! ============================================================================
!                                assemble_F
! ============================================================================
! Call assemble_F: This program assembles all the forces we have calculated
! in Dassemble_2c and Dassemble_3c, and assemble_usr.
          !write (*,*) '  '
          !write (*,*) '  '
          !write (*,*) ' ***************************************************** '
          !write (*,*) ' Assemble all force contributions. '
          call timer_start_i(50)
          call assemble_F ! (natoms, itheory, itheory_xc, igauss, ivdw, iharmonic, ibias, iwrtfpieces)
          call timer_stop_i(50)

! Reassign forces for tolerance testing.
          ftotold = ftotnew
          ftotnew = ftot


! Deallocate Arrays
! ===========================================================================


! Format Statements
! ===========================================================================
100     format (2x, 70('='))


        return
        end subroutine getforces_mcweda

