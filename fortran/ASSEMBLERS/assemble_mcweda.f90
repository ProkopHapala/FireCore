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


! readdata.f90
! Program Description
! ===========================================================================
!       This routine reads the different data file.
!
! ===========================================================================

! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine assemble_mcweda ()
        use debug
        use options
        !use outputs
        !use mpi_main
        use loops
        use neighbor_map
        use configuration
        use interactions
        use energy
        !use md
        use timing
        implicit none

! Argument Declaration and Description
! ===========================================================================
! Output

! Local Parameters and Data Declaration
! ===========================================================================

! Local Variable Declaration and Description
! ===========================================================================
        integer iatom
        integer jatom
        integer ineigh
        integer mbeta
        integer kforce

! Procedure
        if(idebugWrite .gt. 0) write(*,*) "BEGIN assemble_mcweda() itheory", itheory
! ===========================================================================

!        !write (*,*) '  '
!        !write (*,100)
!        !write (*,*) ' Now we are assembling McWeda Hamiltonian. '
!        !write (*,*) '  '

! ===========================================================================
!                             neighbors mapping
! ===========================================================================
! Only do these neighbor calculations if we are on the first step of the
! self-consistent cycle.

! Find all neighbors to atoms in the central cell. An atom is at lattice
! vector l(beta), and basis b(j). We refer to this as (beta,j). An atom in
! the central cell is at (0,j). Find all neighbors to atom (0,i)

! neighn(iatom) = # of neighbors of iatom
! neighj(iatom,ineigh) = j-sub-m, the jatom value of the ineigh'th neighbor.
! neighb(iatom,ineigh) = beta-sub-m, the beta value for the ineigh'th neighbor.
          if (Kscf .eq. 1) then
           call timer_start_i(5)
           if (ifixneigh .eq. 0) then
            call reallocate_neigh ! (nprocs, my_proc, iordern, itheory, itheory_xc, igauss, icluster, ivdw, iwrthampiece, iwrtatom, igrid)
            call neighbors   ! (nprocs, my_proc, iordern, icluster, iwrtneigh, ivdw)
            call neighborsPP ! (nprocs, my_proc, iordern, icluster,iwrtneigh)
            call num_neigh_tot ! (numorb_max)
! bias voltage option
!             if (ibias .eq. 1) call reallocate_bias (natoms)     ! IF_DEF_BIAS_END
           else
            !write (*,*) ' Using neighbor map from NEIGHBORS file. '
            call initneighbors ! (natoms, ivdw, nstepi)
            call num_neigh_tot ! (numorb_max)
           end if
           call backnay ()
            !SFIRE
            call neighbors_pairs ! (icluster)
            !SFIRE
           call common_neighbors   ! (nprocs, my_proc, iordern, iwrtneigh_com)
           call common_neighborsPP ! (nprocs, my_proc, iordern,  iwrtneigh_com, icluster)
           call timer_stop_i(5)

           call check_neighbors()
          end if ! end if (Kscf .eq. 1)
          !call check_neighbors()
! ===========================================================================
!                              ewald energy
! ===========================================================================
          call timer_start_i(6)
          if ((itheory .eq. 1 .or. itheory .eq. 2) .and. Kscf .eq. 1) then
           kforce = 0
           call get_ewald ! (nprocs, my_proc, kforce, icluster, itheory, iordern)
          end if
          call timer_stop_i(6)

! ===========================================================================
! ---------------------------------------------------------------------------
!                   A S S E M B L E    H A M I L T O N I A N
!              A N D   O B T A I N   B A N D - S T R U C T U R E
! ---------------------------------------------------------------------------
! ===========================================================================
! Assemble the matrix elements of the Hamiltonian - all in real space.

! Set up neigh_self.  The variable neigh_self(natoms) is the ineigh value
! for the "self interaction".  Find the neighbor-number of iatom with itself
! (neigh_self) in order to put the result of VNA_atom (doscentros) into
! VNA(mu,nu,iatom,neigh_self).
! Initialize to something ridiculous.

          neigh_self = -999
          do iatom = 1, natoms
           do ineigh = 1, neighn(iatom)
            mbeta = neigh_b(ineigh,iatom)
            jatom = neigh_j(ineigh,iatom)
            if (iatom .eq. jatom .and. mbeta .eq. 0) neigh_self(iatom) = ineigh
           end do
          end do
          neighPP_self = -999
          do iatom = 1, natoms
           do ineigh = 1, neighPPn(iatom)
            mbeta = neighPP_b(ineigh,iatom)
            jatom = neighPP_j(ineigh,iatom)
            if (iatom .eq. jatom .and. mbeta .eq. 0)  neighPP_self(iatom) = ineigh
           end do
          end do

! ! IF_DEF_VDW
!           if (ivdw .eq. 1) then
!            neigh_vdw_self = -999
!            do iatom = 1, natoms
!             do ineigh = 1, neighn_vdw(iatom)
!              mbeta = neigh_b_vdw(ineigh,iatom)
!              jatom = neigh_j_vdw(ineigh,iatom)
!              if (iatom .eq. jatom .and. mbeta .eq. 0)   neigh_vdw_self(iatom) = ineigh
!             end do
!            end do
!           end if
! ! IF_DEF_VDW

! ===========================================================================
!                               assemble_1c
! ===========================================================================
! Assemble the one-center exchange-correlation interactions.
          call timer_start_i(7)
          if(itheory_xc .eq. 2 ) then
           !write (*,*) ' ***************************************************** '
           !write (*,*) ' Assemble one-center interactions. '
           !write (*,*) ' Assemble OLS-xc exchange-correlation interactions. '
           call assemble_olsxc_1c ! (natoms, itheory, iforce)
          endif
         !if (V_intra_dip .eq. 1) call assemble_1c_vdip ! (natoms, itheory, iforce)    ! 
          call timer_stop_i(7)
! ===========================================================================
!                               assemble_2c
! ===========================================================================
! We now do the 2center contributions - assemble_2c does the assembly and
! doscentros does the interpolating and "fundamental" calculations.
! assemble_2c ONLY does 2c terms. No 3c terms allowed. See assemble_3c
! and trescentros for 3c terms.
          call timer_start_i(8)

          if (Kscf .eq. 1) then
           call timer_start_i(20)
           !write (*,*) ' Assemble two-center interactions. '
           call assemble_sVNL()  ! (iforce)
           call assemble_2c()    ! (nprocs, iforce, iordern, ioff2c)
           call assemble_2c_PP() ! (nprocs, iforce, iordern)
! assemble bias matrix
!            if (ibias .eq. 1) call assemble_Vbias ( nprocs, iforce, iordern, ioff2c )   ! IF_DEF_BIAS_END
           call timer_stop_i(20)
          end if ! end if of Kscf = 1

! Call the exchange-correlation interactions based on method chosen
! (i.e. itheory_xc).
          if (itheory_xc .eq. 1 ) then
          !write (*,*) ' Assemble SN-xc exchange-correlation interactions. '
          call timer_start_i(21)
          if (itheory .eq. 1) then
           call average_ca_rho ! (nprocs, Kscf, iforce, iordern, igauss)
          else
           call average_rho ! (nprocs, Kscf, iforce, iordern, igauss)
          endif
          call timer_stop_i(21)

           call timer_start_i(22)
           call assemble_snxc_on ( uxcdcc_sn ) ! (natoms, nprocs, my_proc, iordern, itheory, uxcdcc_sn)
           call timer_stop_i(22)

           call timer_start_i(23)
           call assemble_snxc_off ! (natoms, nprocs, my_proc, iordern,  itheory)
           call timer_stop_i(23)
          end if ! if (itheory_xc = 1)

          if (itheory_xc .eq. 2 ) then
           !write (*,*) ' Assemble OLS-xc exchange-correlation interactions.'
           
           call timer_start_i(21)
           if (itheory .eq. 1) then
            call average_ca_rho ! (nprocs, Kscf, iforce, iordern, igauss)
            !call average_rho (nprocs, Kscf, iforce, iordern, igauss)
           else
            call average_rho ! (nprocs, Kscf, iforce, iordern, igauss)
           endif
           call timer_stop_i(21)

           call timer_start_i(22)
           call assemble_olsxc_on (uxcdcc_ols) ! (natoms, nprocs, my_proc, iordern, itheory, uxcdcc_ols)
           call timer_stop_i(22)

           call timer_start_i(23)
           call assemble_olsxc_off ! (nprocs, my_proc, iordern, itheory)
           call timer_stop_i(23)
          end if ! if (itheory_xc = 2)

!JIMM
          call timer_start_i(24)
          if (itheory .eq. 1) then
           !write (*,*) ' Assemble two-center DOGS interactions. '
           !if (idipole .eq. 1) then                                ! IF_DEF_DIPOL_END
           !     call assemble_ca_2c_dip !(nprocs, iforce, iordern) ! IF_DEF_DIPOL_END
           !else                                                    ! IF_DEF_DIPOL_END
                call assemble_ca_2c     !(nprocs, iforce, iordern)
           ! end if                                                 ! IF_DEF_DIPOL_END
          endif
          call timer_stop_i(24)
          call timer_stop_i(8)

          !call debug_writeBlockedMat( "H_2c_t.log", t_mat )
          !call debug_writeBlockedMat( "H_2c_vna.log", vna )
          !call debug_writeBlockedMat( "H_2c_vxc.log", vxc )
          !call debug_writeBlockedMat( "H_2c_vxc1c.log", vxc_1c )
          !call debug_writeBlockedMat( "H_2c_vca.log", vca )
          !call debug_writeBlockedMat( "H_2c_vxcca.log", vxc_ca )
          !call debug_writeBlockedMat( "H_2c_ewaldLR.log", ewaldlr )
          !call debug_writeBlockedMat( "H_2c_ewaldSR.log", ewaldsr )
          !call debug_writeBlockedMat( "H_mat_2c.log", h_mat )
! ===========================================================================
!                               assemble_3c
! ===========================================================================
! We now do the 3center contributions - assemble_3c does the assembly and
! trecentros does the interpolating and "fundamental" calculations.
! assemble_3c ONLY does 3c terms. No 2c terms allowed. See assemble_2c
! and doscentros for 2c terms.
          
          call timer_start_i(9)
          if (Kscf .eq. 1) then
           
           call timer_start_i(31)
           call assemble_3c    ! (nprocs, iordern, igauss, itheory_xc)
           call timer_stop_i(31)

           call timer_start_i(32)
           call assemble_3c_PP ! (nprocs, iordern)
           call timer_stop_i(32)

! JIMM
! IF_DEF_QMMM
!           if (iqmmm .eq.1 ) then
!             !write (*,*) ' Assemble qm/mm interactions. '
!             if (idipole .eq. 0) call assemble_qmmm (nprocs, iordern)
!             if (idipole .eq. 1) call assemble_qmmm_dip (nprocs, iordern)
!           else
!             eqmmm = 0.0d0
!             ewaldqmmm = 0.0d0
!           end if
! END_DEF_QMMM
!JIMM

          end if

          if (itheory .eq. 1) then
           call timer_start_i(33)
           !write (*,*) ' Assemble three-center DOGS interactions. '
           ! if (idipole .eq. 1) then                              ! IF_DEF_DIPOL_END
           !     call assemble_ca_3c_dip (nprocs, iordern, igauss) ! IF_DEF_DIPOL_END
           ! else                                                  ! IF_DEF_DIPOL_END
                call assemble_ca_3c ! (nprocs, iordern, igauss)
           !end if                                                 ! IF_DEF_DIPOL_END
           call timer_stop_i(33)

           call timer_start_i(34)
! Add assemble_lr here for the long long-range ewald contributions
           !write (*,*) ' Assemble long-range interactions. '
           !if (idipole .eq. 1)                              ! IF_DEF_DIPOL_END
           !     call assemble_lr_dip (nprocs, iordern)      ! IF_DEF_DIPOL_END
           !else                                             ! IF_DEF_DIPOL_END
                call assemble_lr ! (nprocs, iordern)
           !end if                                           ! IF_DEF_DIPOL_END
           call timer_stop_i(34)     
          endif
          call timer_stop_i(9)
          !call debug_writeBlockedMat( "H_3c_t.log", t_mat )
          !call debug_writeBlockedMat( "H_3c_vna.log", vna )
          !call debug_writeBlockedMat( "H_3c_vxc.log", vxc )
          !call debug_writeBlockedMat( "H_2c_vxc1c.log", vxc_1c )
          !call debug_writeBlockedMat( "H_3c_vca.log", vca )
          !call debug_writeBlockedMat( "H_3c_vxcca.log", vxc_ca )
          !call debug_writeBlockedMat( "H_3c_ewaldLR.log", ewaldlr )
          !call debug_writeBlockedMat( "H_3c_ewaldSR.log", ewaldsr )
          !call debug_writeBlockedMat( "H_mat_3c.log", h_mat )

          !write (*,*) ' ***************************************************** '

! ===========================================================================
!                                 Build H
! ===========================================================================
! Set up the full Hamiltonian and !writeout HS.dat.
          call timer_start_i(10)
          call buildh !(nprocs, itheory, iordern, itestrange, testrange, ibias, iwrtHS)
          call timer_stop_i(10)
! ===========================================================================
! For iwrthampiece .eq. 1 (file - output.input), !write out Hamiltonian pieces
!          if (iwrthampiece .eq. 1) call hampiece (itheory)  ! IF_DEF_HAMPIECES_END


! Deallocate Arrays
! ===========================================================================

! Format Statements
! ===========================================================================
100     format (2x, 70('='))


        return
        end subroutine assemble_mcweda

