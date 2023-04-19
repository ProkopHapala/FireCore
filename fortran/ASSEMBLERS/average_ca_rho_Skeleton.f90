      subroutine average_ca_rho () ! (nprocs, Kscf, iforce, iordern, igauss)
!*****************************************************************************
!  1.    ON-SITE PART 
!*****************************************************************************
         do iatom = iatomstart, iatomstart - 1 + natomsp
            call doscentrosS (interaction0, isorp, iforce, in1, in2, in3, y, eps, sm, spm)
            do ineigh = 1, neighn(iatom)
               call epsilon   (r2, sighat, eps)
               call deps2cent (r1, r2, eps, deps)
               if (iatom .eq. jatom .and. mbeta .eq. 0) then
                  do isorp = 1, nssh(in2)
                     call doscentros  (17, isorp, iforce, in1, in2, in1, y, eps, deps, rhomx, rhompx)
                     call doscentrosS (22, isorp, iforce, in1, in2, in1, y, eps,       rhomm, rhompm)
                     ! UPDATE_RHO()
                  end do ! endo do isorp
               else
                  do isorp = 1, nssh(in2)
                     call doscentros (17, isorp, iforce, in1, in2,  in1, y, eps, deps, rhomx, rhompx)
                     call doscentrosS(22, isorp, iforce, in1, in2,  in1, y, eps,       rhomm, rhompm)
                     ! UPDATE_RHO()
                  end do ! endo do isorp
                  ! UPDATE_RHO()
               end if ! end if (iatom.eq.jatom)
            end do !!!!! end do ineigh
            ! UPDATE_RHO()
         enddo !!!!! end do iatom    

        
!*****************************************************************************
!  2.    OFF-SITE PART
!  2.1.     THREE-CENTER PART  of   OFF-SITE PART 
!*****************************************************************************
         do ialp = iatomstart, iatomstart - 1 + natomsp   ! Loop over the atoms in the central cell.
            do ineigh = 1, neigh_comn(ialp)
               if (mneigh .ne. 0) then
                  call epsilon (rhat, sighat, eps)
                  do isorp = 1, nssh(indna)
                     call trescentros (3, isorp, isorpmax, in1, in2, indna, x, y, cost, eps, rhomx, nspecies)
                     call trescentrosS(   isorp, isorpmax, in1, in2, indna, x, y, cost, eps, rhomm, nspecies)
                     ! UPDATE_RHO()
                  end do ! do isorp
               end if
            end do ! end do neigh_comn(ialp)
         end do  ! end do ialp


! ****************************************************************************
!  2.2.  TWO-CENTER PART
! ****************************************************************************
         do iatom = iatomstart, iatomstart - 1 + natomsp
            do ineigh = 1, neighn(iatom)          ! <==== loop over i's neighbors
               if (iatom .eq. jatom .and. mbeta .eq. 0) then
               ! Do nothing here - special case. Interaction already calculated in on-site stuff, i.e. assemble_olsxc_on.f90 We calculate only off diagonal elements here.
               else
                  call epsilon (r2, sighat, eps)
                  call deps2cent (r1, r2, eps, deps)
                  ! A. Left piece: den_ontopl <i|n_i|j> (den1,2) part of <i|n|j> (denmx)
                  do isorp = 1, nssh(in3)
                     call doscentros (15, isorp, iforce, in1, in3, in2, y, eps, deps, rhomx, rhompx)
                     call doscentrosS(20, isorp, iforce, in1, in3, in2, y, eps,       rhomm, rhompm)
                     ! UPDATE_RHO()
                  end do
                  ! B. Right piece: den_ontopr <i|n_j|j> (den0) part of <i|n|j> (denmx)
                  do isorp = 1, nssh(in3)
                     call doscentros (16, isorp, iforce, in1, in2, in2, y, eps, deps, rhomx, rhompx)
                     call doscentrosS(21, isorp, iforce, in1, in2, in2, y, eps,       rhomm, rhompm )
                     ! UPDATE_RHO()
                  end do
                  if (Kscf .eq. 1) then
                     call doscentrosS (23, 0, iforce, in1, in2, in2, y, eps, sm, spm )
                     ! UPDATE_S()
                  else
                    ! UPDATE_S()
                  end if
                  ! UPDATE_RHO()
               end if
            end do ! do ineigh
         end do  ! do iatom

! Deallocate Arrays
! ===========================================================================
         deallocate (rhom_3c)

! Format Statements
! ===========================================================================
  100    format(9f14.6)
         return
      end subroutine average_ca_rho
