! copyright info:
!
!                             @Copyright 2002
!                           Fireball Committee
! Brigham Young University - James P. Lewis, Chair
! Arizona State University - Otto F. Sankey
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


! buildh.f90
! Program Description
! ===========================================================================
!       This routines combines all of the interations into one Hamiltonian
! matirx.
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
        subroutine buildh ! ( nprocs, itheory, iordern, itestrange, testrange, ibias, iwrtHS )
        use options
        use configuration
        use interactions
        use neighbor_map
        use dimensions
        use density
        ! use bias
        !use options, only: V_intra_dip
        implicit none

! Argument Declaration and Description
! ===========================================================================
        ! integer, intent (in) :: iordern
        ! integer, intent (in) :: itestrange
        ! integer, intent (in) :: itheory
        ! integer, intent (in) :: nprocs
        ! integer, intent (in) :: ibias
        ! integer, intent (in) :: iwrtHS
        ! real,    intent (in) :: testrange

! Local Parameters and Data Declaration
! ===========================================================================

! Local Variable Declaration and Description
! ===========================================================================
        integer katom
        integer iatom
        integer iatomstart
        integer ierror
        integer imu
        integer in1
        integer in2
        integer ineigh
        integer matom
        integer inu
        integer jatom
        integer mbeta
        integer my_proc
        integer natomsp

        integer nnu, nmu

        real distance

        real, dimension (numorb_max, numorb_max) :: htemp
        real, dimension (numorb_max, numorb_max) :: stemp

        real, dimension (3) :: dvec

        integer issh
        integer numorb
        integer jatom0
        integer ineigh0
        integer mbeta0

! BTN communication domain
!        integer MPI_BTN_WORLD, MPI_OPT_WORLD, MPI_BTN_WORLD_SAVE                  ! IF_DEF_ORDERN_END
!        common  /btnmpi/ MPI_BTN_WORLD, MPI_OPT_WORLD, MPI_BTN_WORLD_SAVE         ! IF_DEF_ORDERN_END

! Allocate Arrays
! ===========================================================================

! Procedure
! ===========================================================================
! Find out which processor this is.
! ! IF_DEF_ORDERN
!        if (iordern .eq. 1) then
!         call MPI_COMM_RANK (MPI_BTN_WORLD, my_proc, ierror)
!         natomsp = natoms/nprocs
!         if (my_proc .lt. mod(natoms,nprocs)) then
!          natomsp = natomsp + 1
!          iatomstart = natomsp*my_proc + 1
!         else
!          iatomstart = (natomsp + 1)*mod(natoms,nprocs)  + natomsp*(my_proc - mod(natoms,nprocs)) + 1
!         end if
!        else
! ! END_DEF_ORDERN
         iatomstart = 1
         natomsp = natoms
!        end if  ! IF_DEF_ORDERN_END


! Set up the full Hamiltonian.
!$omp parallel do private (in1, ineigh, mbeta, jatom, in2, inu, imu, distance)
!!$omp parallel do private (in1, in2, mbeta, jatom, distance) --this is the old setting
        do iatom = iatomstart, iatomstart - 1 + natomsp
         in1 = imass(iatom)
         do ineigh = 1, neighn(iatom)
          mbeta = neigh_b(ineigh,iatom)
          jatom = neigh_j(ineigh,iatom)
          in2 = imass(jatom)

          nnu = num_orb(in2)
          nmu = num_orb(in1)

! Here we add in all of the charge interactions for Kohn-Sham method .
          if (itheory .eq. 3) then
           !do inu = 1, num_orb(in2)
            !do imu = 1, num_orb(in1)
              !h_mat(imu,inu,ineigh,iatom) =  t_mat(imu,inu,ineigh,iatom) + vna(imu,inu,ineigh,iatom) + vxc(imu,inu,ineigh,iatom) + vca(imu,inu,ineigh,iatom)
            !end do
           !end do

        h_mat(:nmu,:nnu,ineigh,iatom) =  t_mat(:nmu,:nnu,ineigh,iatom) + vna(:nmu,:nnu,ineigh,iatom) + vxc(:nmu,:nnu,ineigh,iatom) + vca(:nmu,:nnu,ineigh,iatom)

          else

! Add all the pieces together for the Hamiltonian.
           !do inu = 1, num_orb(in2)
            !do imu = 1, num_orb(in1)
             !h_mat(imu,inu,ineigh,iatom) = t_mat(imu,inu,ineigh,iatom) + vna(imu,inu,ineigh,iatom) + vxc(imu,inu,ineigh,iatom) + vxc_1c(imu,inu,ineigh,iatom)
            !end do
           !end do
        
        h_mat(:nmu,:nnu,ineigh,iatom) = t_mat(:nmu,:nnu,ineigh,iatom) + vna(:nmu,:nnu,ineigh,iatom) + vxc(:nmu,:nnu,ineigh,iatom) + vxc_1c(:nmu,:nnu,ineigh,iatom)

! Here we add in all of the charge interactions, as well as the long range
! coulomb(ewaldlr) and the short range correction, 1/R (ewaldsr).
           if (itheory .eq. 1 .or. itheory .eq. 2) then
            !do inu = 1, num_orb(in2)
             !do imu = 1, num_orb(in1)
!$omp atomic
              !h_mat(imu,inu,ineigh,iatom) = h_mat(imu,inu,ineigh,iatom) + vca(imu,inu,ineigh,iatom) + vxc_ca(imu,inu,ineigh,iatom) + ewaldlr(imu,inu,ineigh,iatom) - ewaldsr(imu,inu,ineigh,iatom)
            !   h_mat(imu,inu,ineigh,iatom) = h_mat(imu,inu,ineigh,iatom) + ewaldqmmm(imu,inu,ineigh,iatom)   ! IF_DEF_QMMM_END
             !end do ! do imu
            !end do ! do inu

            h_mat(:nmu,:nnu,ineigh,iatom) = h_mat(:nmu,:nnu,ineigh,iatom) + vca(:nmu,:nnu,ineigh,iatom) + vxc_ca(:nmu,:nnu,ineigh,iatom) + ewaldlr(:nmu,:nnu,ineigh,iatom) - ewaldsr(:nmu,:nnu,ineigh,iatom)

           end if ! if (itheory .eq. 1)
          endif ! if (iKS .eq. 1)
         end do ! do ineigh

! ! IF_DEF_IXCZW
!             !New V_intra_dip_1c JUNE 2019
!             if (V_intra_dip .eq. 1) then
!               matom = neigh_self(iatom)
!               do inu = 1, num_orb(in1)
!                do imu = 1, num_orb(in1)
!                 h_mat(imu,inu,matom,iatom) = h_mat(imu,inu,matom,iatom) + Vdip_1c(imu,inu,iatom)
!                end do ! do imu
!               end do ! do inu
!             end if ! if (V_intra_dip .eq. 1) 
!             !End of New V_intra_dip_1c JUNE 2019
! ! END_DEF_IXCZW

        end do ! do iatom

! ! IF_DEF_BIAS
! ! Here we add extra terms of Hamiltonian comming from the bias voltage (not KS)
!         if (ibias .eq. 1 .and. itheory .ne. 3) then
!          do iatom = 1, natoms
!           in1 = imass(iatom)
! ! loop over neighbors
!           do ineigh = 1, neighn(iatom)
!            mbeta = neigh_b(ineigh,iatom)
!            jatom = neigh_j(ineigh,iatom)
!            in2 = imass(jatom)
! ! Add all the pieces together for the Hamiltonian.
!            do inu = 1, num_orb(in2)
!             do imu = 1, num_orb(in1)
!              h_mat(imu,inu,ineigh,iatom) = h_mat(imu,inu,ineigh,iatom)       + Vbias_mat(imu,inu,ineigh,iatom)
!             end do
!            end do
!           enddo ! do ineigh
!          enddo ! do iatom
!         endif ! if(ibias)
! ! END_DEF_BIAS


! We set matrix elements equal to zero if outside our test range.
        do iatom = iatomstart, iatomstart - 1 + natomsp
         in1 = imass(iatom)
         do ineigh = 1, neighn(iatom)
          mbeta = neigh_b(ineigh,iatom)
          jatom = neigh_j(ineigh,iatom)
          in2 = imass(jatom)
          if (itestrange .eq. 0) then
           distance =                                                        &
     &      sqrt((ratom(1,iatom) - (xl(1,mbeta) + ratom(1,jatom)))**2        &
     &           + (ratom(2,iatom) - (xl(2,mbeta) + ratom(2,jatom)))**2      &
     &           + (ratom(3,iatom) - (xl(3,mbeta) + ratom(3,jatom)))**2)
           if (distance .gt. testrange) then
                nnu = num_orb(in2)
                nmu = num_orb(in1)
                
                
            !do inu = 1, num_orb(in2)
             !do imu = 1, num_orb(in1)
              !h_mat(imu,inu,ineigh,iatom) = 0.0d0
              !t_mat(imu,inu,ineigh,iatom) = 0.0d0
              !s_mat(imu,inu,ineigh,iatom) = 0.0d0
              !vna(imu,inu,ineigh,iatom) = 0.0d0
              !vxc(imu,inu,ineigh,iatom) = 0.0d0
              !if (itheory .eq. 1 .or. itheory .eq. 2) then
              ! ewaldlr(imu,inu,ineigh,iatom) = 0.0d0
              ! ewaldsr(imu,inu,ineigh,iatom) = 0.0d0
!             !  ewaldqmmm(imu,inu,ineigh,iatom) = 0.0d0    ! IF_DEF_QMMM_END
              !end if
             !end do
            !end do

              h_mat(:nmu,:nnu,ineigh,iatom) = 0.0d0
              t_mat(:nmu,:nnu,ineigh,iatom) = 0.0d0
              s_mat(:nmu,:nnu,ineigh,iatom) = 0.0d0
              vna  (:nmu,:nnu,ineigh,iatom) = 0.0d0
              vxc  (:nmu,:nnu,ineigh,iatom) = 0.0d0
              if (itheory .eq. 1 .or. itheory .eq. 2) then
                !ewaldlr(:nmu,:nnu,ineigh,iatom) = 0.0d0
                ewaldsr(:nmu,:nnu,ineigh,iatom) = 0.0d0
!               !ewaldqmmm(:nmu,:nnu,ineigh,iatom) = 0.0d0    ! IF_DEF_QMMM_END
              end if

           end if
          end if
         end do ! do ineigh
        end do ! do iatom


!        if (iordern .eq. 1) call buildh_ordern_final (natoms, nprocs, my_proc, itheory)    ! IF_DEF_ORDERN_END

! ===========================================================================
!         if (iwrtHS .eq. 1) then
! ! Writeout HS.dat
!         write (*,*)
!         write (*,*) ' The information needed to make a TB model or '
!         write (*,*) ' for doing the complex band-structure is in '
!         write (*,*) ' the file ---> HS.dat'

!         open (unit = 11, file = 'HS.dat', status = 'unknown')

!         numorb = 0
!         do iatom = 1, natoms
!          in1 = imass(iatom)
!          do ineigh = 1, neighn_tot(iatom)
!           jatom = neighj_tot(ineigh,iatom)
!           in2 = imass(jatom)
!           do imu = 1, num_orb(in1)
!            do inu = 1, num_orb(in2)
!               numorb = numorb + 1
!            end do  ! do imu
!           end do  ! do inu
! !
!          end do  ! do ineigh
!         end do ! do iatom
! !
!         write (11,*) numorb

!         do iatom = 1, natoms
!          in1 = imass(iatom)
! ! loop over total list of neighbors
!          do ineigh = 1, neighn_tot(iatom)
!           jatom = neighj_tot(ineigh,iatom)
!           mbeta = neighb_tot(ineigh,iatom)
!           in2 = imass(jatom)
!           dvec(:) = ratom(:,jatom) + xl(:,mbeta) - ratom(:,iatom)
!           htemp = 0.0d0
!           stemp = 0.0d0

! ! loop over regular list of neighbors
!           do ineigh0 = 1, neighn(iatom)
!            jatom0 = neigh_j(ineigh0,iatom)
!            mbeta0 = neigh_b(ineigh0,iatom)

! ! find identical neighbors
!            if (jatom .eq. jatom0 .and. mbeta .eq. mbeta0) then
!             do imu = 1, num_orb(in1)
!              do inu = 1, num_orb(in2)
!               htemp(imu,inu) = htemp(imu,inu) + h_mat(imu,inu,ineigh0,iatom)
!               stemp(imu,inu) = s_mat(imu,inu,ineigh0,iatom)
!              enddo ! inu
!             enddo ! imu
!            endif
!           enddo ! do ineigh0
! ! loop over PP list of neighbors
!           do ineigh0 = 1, neighPPn(iatom)
!            jatom0 = neighPP_j(ineigh0,iatom)
!            mbeta0 = neighPP_b(ineigh0,iatom)

! ! find identical neighbors
!            if (jatom .eq. jatom0 .and. mbeta .eq. mbeta0) then
!             do imu = 1, num_orb(in1)
!              do inu = 1, num_orb(in2)
!               htemp(imu,inu) = htemp(imu,inu) + vnl(imu,inu,ineigh0,iatom)
!              enddo ! inu
!             enddo ! imu
!            endif
!           enddo ! do ineigh0

! ! write elements
!           do imu = 1, num_orb(in1)
!            do inu = 1, num_orb(in2)
!             write (11,50) iatom, imu, jatom, inu, htemp(imu,inu), stemp(imu,inu), dvec
!            end do  ! do imu
!           end do  ! do inu
! !
!          end do  ! do ineigh
!         end do ! do iatom

!         close (unit = 11)
!         end if ! (iwrtHS .eq. 1)

! Format Statements
! ===========================================================================
 50     format (4i4, 2f12.6, 3f12.6)

        return
      end subroutine buildh
