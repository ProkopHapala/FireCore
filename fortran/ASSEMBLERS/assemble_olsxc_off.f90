! copyright info:
!
!                             @Copyright 2005
!                           Fireball Committee
! Brigham Young University - James P. Lewis, Chair
! Arizona State University - Otto F. Sankey
! Universidad de Madrid - Jose Ortega
! Universidad de Madrid - Pavel Jelinek

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

 
! assemble_olsxc_2c.f90
! Program Description
! ===========================================================================
!       This routine assembles the two & three-center exchange-correlation
! for the average density approximation. 
!
! ===========================================================================
! Code written by:
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
        subroutine assemble_olsxc_off ! (nprocs, my_proc, iordern, itheory)
        use options
        use configuration
        use charges
        use density
        use dimensions
        use debug
        use interactions
        use neighbor_map
        use timing
        implicit none
 
! Argument Declaration and Description
! ===========================================================================
! Input
        ! integer, intent (in) :: iordern
        ! integer, intent (in) :: itheory
        ! integer, intent (in) :: my_proc
        ! integer, intent (in) :: nprocs


! Local Parameters and Data Declaration
! ===========================================================================
 
! Local Variable Declaration and Description
! ===========================================================================
        integer iatom
        integer iatomstart
        integer imu
        integer in1, in2, in3
        integer ineigh
        integer interaction
        integer inu
        integer isorp
        integer jatom
        integer kforce
        integer matom
        integer mbeta
        integer natomsp

        real  y 
        real dxn
        real, dimension (numorb_max, numorb_max) :: bcxcx
        real, dimension (numorb_max, numorb_max) :: denmx
        real, dimension (numorb_max, numorb_max) :: den1x
        real, dimension (numorb_max, numorb_max) :: rhomx
        real, dimension (3, numorb_max, numorb_max) :: rhompx
        real, dimension (3, 3) :: eps
        real, dimension (3, 3, 3) :: deps
        real, dimension (3) :: r1, r2, r21
        real, dimension (3) :: sighat
        real, dimension (numorb_max, numorb_max) :: sx
        real exc, muxc, dexc, d2exc, dmuxc, d2muxc
        real excij, muxcij, dexcij, d2excij, dmuxcij, d2muxcij
        integer issh, jssh

! Procedure
! ===========================================================================
! Determine which atoms are assigned to this processor.
! ! IF_DEF_ORDERN
!         if (iordern .eq. 1) then
!          natomsp = natoms/nprocs
!          if (my_proc .lt. mod(natoms,nprocs)) then
!           natomsp = natomsp + 1
!           iatomstart = natomsp*my_proc + 1
!          else
!           iatomstart = (natomsp + 1)*mod(natoms,nprocs)                      &
!            + natomsp*(my_proc - mod(natoms,nprocs)) + 1
!          end if
!         else
! ! END_DEF_ORDERN
         iatomstart = 1
         natomsp = natoms
!        end if  ! IF_DEF_ORDERN_END

         ncall_assemble_olsxc_off=ncall_assemble_olsxc_off+1

! Loop over the atoms in the central cell.
!$omp parallel do private ()
        do iatom = iatomstart, iatomstart - 1 + natomsp
         matom = neigh_self(iatom)
         r1(:) = ratom(:,iatom)
         in1 = imass(iatom)

! Loop over the neighbors of each iatom.
         do ineigh = 1, neighn(iatom)       ! <==== loop over i's neighbors
          mbeta = neigh_b(ineigh,iatom)
          jatom = neigh_j(ineigh,iatom)
          r2(:) = ratom(:,jatom) + xl(:,mbeta)
          in2 = imass(jatom)

! SET-UP STUFF
! ****************************************************************************
! Find r21 = vector pointing from r1 to r2, the two ends of the bondcharge.
! This gives us the distance dbc (or y value in the 2D grid).
          r21(:) = r2(:) - r1(:)
          y = sqrt(r21(1)*r21(1) + r21(2)*r21(2) + r21(3)*r21(3))
 
! Find the unit vector in sigma direction.
          if (y .lt. 1.0d-05) then
           sighat(1) = 0.0d0
           sighat(2) = 0.0d0
           sighat(3) = 1.0d0
          else
           sighat(:) = r21(:)/y
          end if
              
          call epsilon (r2, sighat, eps)
          call deps2cent (r1, r2, eps, deps)

! --------------------------
! DEBUG : TO EXPORT For checking /pyBall/FireballOCL/OCL_Hamiltonian.py
! DEBUG OCL PARITY START (assemble_olsxc_off)
          if (idebugWrite .gt. 0 .and. iatom .eq. 2 .and. ineigh .eq. 1) then
             write(*,*) '[XC_OFF][geom] iatom,ineigh,jatom,mbeta,in1,in2,y=', iatom, ineigh, jatom, mbeta, in1, in2, y
          end if
! DEBUG OCL PARITY END
! -------------------------- 
 
! ****************************************************************************
!
! CALL DOSCENTROS AND GET VXC FOR ATM CASE - AVERAGE DENSITY APPROXIMATION 
! ****************************************************************************
! If r1 .ne. r2, then this is a case where the potential is located
! at one of the sites of a wavefunction (ontop case).
          if (iatom .eq. jatom .and. mbeta .eq. 0) then
 
! Do nothing here - special case. Interaction already calculated in atm case.
 
          else

! This is the ontop case for the exchange-correlation energy
! Horsfield like term <i mu|Vxc(rho_i+j)| j nu>
           isorp = 0
           interaction = 6
           in3 = in2
           call doscentros (interaction, isorp, kforce, in1, in2, in3, y,     eps, deps, rhomx, rhompx)
! --------------------------
! DEBUG : TO EXPORT For checking /pyBall/FireballOCL/OCL_Hamiltonian.py           
! DEBUG OCL PARITY START (assemble_olsxc_off)
           if (idebugWrite .gt. 0 .and. iatom .eq. 2 .and. ineigh .eq. 1) then
              write(*,*) '[XC_OFF][doscentros] interaction=6 iatom,ineigh=', iatom, ineigh
              write(*,'(a,4(1x,e16.8))') '  rhomx row1:', rhomx(1,1), rhomx(1,2), rhomx(1,3), rhomx(1,4)
              write(*,'(a,4(1x,e16.8))') '  rhomx row2:', rhomx(2,1), rhomx(2,2), rhomx(2,3), rhomx(2,4)
              write(*,'(a,4(1x,e16.8))') '  rhomx row3:', rhomx(3,1), rhomx(3,2), rhomx(3,3), rhomx(3,4)
              write(*,'(a,4(1x,e16.8))') '  rhomx row4:', rhomx(4,1), rhomx(4,2), rhomx(4,3), rhomx(4,4)
           end if
! DEBUG OCL PARITY END
! --------------------------           
           do inu = 1, num_orb(in3)
            do imu = 1, num_orb(in1)
!$omp atomic
             vxc(imu,inu,ineigh,iatom) =    vxc(imu,inu,ineigh,iatom) + rhomx(imu,inu)
            end do
           end do

!dani.JOM[
          if ( itheory.eq.1 ) then
! the vxc_ontopl case, transfer of charge (McWeda)
           interaction = 18
           in3 = in2
           do isorp = 1, nssh(in1)
            call doscentros (interaction, isorp, kforce, in1, in1, in3, y, eps, deps, rhomx, rhompx)
            dxn = (Qin(isorp,iatom) - Qneutral(isorp,in1))
            do inu = 1, num_orb(in3)
             do imu = 1, num_orb(in1)
             vxc_ca(imu,inu,ineigh,iatom) =  vxc_ca(imu,inu,ineigh,iatom) + rhomx(imu,inu)*dxn
             end do
            end do
           end do


! the vxc_ontopr case, transfer of charge (McWeda)
           interaction = 19
           in3 = in2
           do isorp = 1, nssh(in2)
            call doscentros (interaction, isorp, kforce, in1, in2, in3, y,  eps, deps, rhomx, rhompx)
            dxn = (Qin(isorp,jatom) - Qneutral(isorp,in2))
            do inu = 1, num_orb(in3)
             do imu = 1, num_orb(in1)
              vxc_ca(imu,inu,ineigh,iatom) = vxc_ca(imu,inu,ineigh,iatom) + rhomx(imu,inu)*dxn
             end do
            end do
           end do
          endif
!dani.JOM]
!

! Restore density and overlap matrices
!  den1x    .... <i|n_i|j> + <i|n_j|j> (on-top interactions)
!  denmx    .... <i|n|j> = <i|n_i|j> + <i|n_j|j> + S_{i.ne.j.ne.k}<i|n_k|j>
!  sx       .... <i|j>
           do inu = 1, num_orb(in3)
            do imu = 1, num_orb(in1)
!$omp atomic
             denmx(imu,inu) = rho_off(imu,inu,ineigh,iatom)
             den1x(imu,inu) = rhoij_off(imu,inu,ineigh,iatom)
             sx(imu,inu) = s_mat(imu,inu,ineigh,iatom)
            end do
           end do

! --------------------------
! DEBUG : TO EXPORT For checking /pyBall/FireballOCL/OCL_Hamiltonian.py  
! DEBUG OCL PARITY START (assemble_olsxc_off)
           if (idebugWrite .gt. 0 .and. iatom .eq. 2 .and. ineigh .eq. 1) then
              write(*,*) '[XC_OFF][inputs] iatom,ineigh= ', iatom, ineigh
              write(*,'(a,4(1x,e16.8))') '  denmx row1:', denmx(1,1), denmx(1,2), denmx(1,3), denmx(1,4)
              write(*,'(a,4(1x,e16.8))') '  denmx row2:', denmx(2,1), denmx(2,2), denmx(2,3), denmx(2,4)
              write(*,'(a,4(1x,e16.8))') '  denmx row3:', denmx(3,1), denmx(3,2), denmx(3,3), denmx(3,4)
              write(*,'(a,4(1x,e16.8))') '  denmx row4:', denmx(4,1), denmx(4,2), denmx(4,3), denmx(4,4)
              write(*,'(a,4(1x,e16.8))') '  den1x row1:', den1x(1,1), den1x(1,2), den1x(1,3), den1x(1,4)
              write(*,'(a,4(1x,e16.8))') '  den1x row2:', den1x(2,1), den1x(2,2), den1x(2,3), den1x(2,4)
              write(*,'(a,4(1x,e16.8))') '  den1x row3:', den1x(3,1), den1x(3,2), den1x(3,3), den1x(3,4)
              write(*,'(a,4(1x,e16.8))') '  den1x row4:', den1x(4,1), den1x(4,2), den1x(4,3), den1x(4,4)
              write(*,'(a,4(1x,e16.8))') '  sx row1:', sx(1,1), sx(1,2), sx(1,3), sx(1,4)
              write(*,'(a,4(1x,e16.8))') '  sx row2:', sx(2,1), sx(2,2), sx(2,3), sx(2,4)
              write(*,'(a,4(1x,e16.8))') '  sx row3:', sx(3,1), sx(3,2), sx(3,3), sx(3,4)
              write(*,'(a,4(1x,e16.8))') '  sx row4:', sx(4,1), sx(4,2), sx(4,3), sx(4,4)

            ! Populate debug buffers
              dbg_vxc_denmx = denmx(1:4,1:4)
              dbg_vxc_den1x = den1x(1:4,1:4)
              dbg_vxc_sx = sx(1:4,1:4)
              ! Store bcxcx after build_olsxc_off call
              call build_olsxc_off (in1, in2, den1x, denmx, sx, ineigh, iatom, bcxcx)
              dbg_vxc_bcxcx = bcxcx(1:4,1:4)
           end if
! DEBUG OCL PARITY END
! --------------------------
           
! Calculate <i| V_xc(n) |j> and <i|V_xc(n_i+n_j)|j>              
           call build_olsxc_off (in1, in2, den1x, denmx, sx, ineigh, iatom,  bcxcx)

! now complete 'non-diagonal' terms <i|V(n)|j>
           if (itheory .eq. 0) then
            do inu = 1, num_orb(in2)
             do imu = 1, num_orb(in1)
              vxc(imu,inu,ineigh,iatom) = vxc(imu,inu,ineigh,iatom) + bcxcx(imu,inu)
             end do
            end do
           else
            do inu = 1, num_orb(in2)
             do imu = 1, num_orb(in1)
              vxc_ca(imu,inu,ineigh,iatom) = vxc_ca(imu,inu,ineigh,iatom) + bcxcx(imu,inu)
             end do
            end do
           end if

! --------------------------
! DEBUG : TO EXPORT For checking /pyBall/FireballOCL/OCL_Hamiltonian.py      
! DEBUG OCL PARITY START (assemble_olsxc_off)
           if (idebugWrite .gt. 0 .and. iatom .eq. 2 .and. ineigh .eq. 1) then
              if (itheory .eq. 0) then
                 write(*,*) '[XC_OFF][vxc] final block iatom,ineigh=', iatom, ineigh
                 write(*,'(a,4(1x,e16.8))') '  vxc row1:', vxc(1,1,ineigh,iatom), vxc(1,2,ineigh,iatom), vxc(1,3,ineigh,iatom), vxc(1,4,ineigh,iatom)
                 write(*,'(a,4(1x,e16.8))') '  vxc row2:', vxc(2,1,ineigh,iatom), vxc(2,2,ineigh,iatom), vxc(2,3,ineigh,iatom), vxc(2,4,ineigh,iatom)
                 write(*,'(a,4(1x,e16.8))') '  vxc row3:', vxc(3,1,ineigh,iatom), vxc(3,2,ineigh,iatom), vxc(3,3,ineigh,iatom), vxc(3,4,ineigh,iatom)
                 write(*,'(a,4(1x,e16.8))') '  vxc row4:', vxc(4,1,ineigh,iatom), vxc(4,2,ineigh,iatom), vxc(4,3,ineigh,iatom), vxc(4,4,ineigh,iatom)
              else
                 write(*,*) '[XC_OFF][vxc_ca] final block iatom,ineigh=', iatom, ineigh
                 write(*,'(a,4(1x,e16.8))') '  vxc_ca row1:', vxc_ca(1,1,ineigh,iatom), vxc_ca(1,2,ineigh,iatom), vxc_ca(1,3,ineigh,iatom), vxc_ca(1,4,ineigh,iatom)
                 write(*,'(a,4(1x,e16.8))') '  vxc_ca row2:', vxc_ca(2,1,ineigh,iatom), vxc_ca(2,2,ineigh,iatom), vxc_ca(2,3,ineigh,iatom), vxc_ca(2,4,ineigh,iatom)
                 write(*,'(a,4(1x,e16.8))') '  vxc_ca row3:', vxc_ca(3,1,ineigh,iatom), vxc_ca(3,2,ineigh,iatom), vxc_ca(3,3,ineigh,iatom), vxc_ca(3,4,ineigh,iatom)
                 write(*,'(a,4(1x,e16.8))') '  vxc_ca row4:', vxc_ca(4,1,ineigh,iatom), vxc_ca(4,2,ineigh,iatom), vxc_ca(4,3,ineigh,iatom), vxc_ca(4,4,ineigh,iatom)
              end if
           end if
! DEBUG OCL PARITY END
! --------------------------

! End if for r1 .ne. r2 case
          end if
! ****************************************************************************
! End loop over iatom and its neighbors - jatom.
         end do
        end do

! Format Statements
! ===========================================================================

        return
        end subroutine assemble_olsxc_off
