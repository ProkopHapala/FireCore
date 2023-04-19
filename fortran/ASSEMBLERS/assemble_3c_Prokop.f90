
        subroutine assemble_3c ! (nprocs, iordern, igauss, itheory_xc)
        use options
        use configuration
        use constants_fireball
        use dimensions
        use forces
        use interactions
        use neighbor_map
        implicit none
 
! Argument Declaration and Description
! ===========================================================================
! Input
        ! integer, intent (in) :: itheory_xc
        ! integer, intent (in) :: igauss
        ! integer, intent (in) :: iordern
        ! integer, intent (in) :: nprocs
 
 
! Local Parameters and Data Declaration
! ===========================================================================
 
! Local Variable Declaration and Description
! ===========================================================================
        integer ialp
        integer iatom
        integer iatomstart
        integer ibeta
        integer ierror
        integer imu
        integer in1
        integer in2
        integer indna
        integer ineigh
        integer interaction
        integer inu
        integer isorp
        integer issh
        integer jatom
        integer jbeta
        integer jssh
        integer mneigh
        integer my_proc
        integer natomsp
        integer jneigh       

        real cost
        real distance_13
        real distance_23
        real dstn_temp
        real stn_temp1
        real stn_temp2
        real x
        real y

        real, dimension (numorb_max, numorb_max) :: bcnax
        real, dimension (3, 3) :: eps
        real, dimension (3) :: r1
        real, dimension (3) :: r2
        real, dimension (3) :: r21
        real, dimension (3) :: r31
        real, dimension (3) :: r32
        real, dimension (3) :: r13
        real, dimension (3) :: r23
        real, dimension (3) :: rhat
        real, dimension (3) :: rna
        real, dimension (3) :: rnabc
        real, dimension (3) :: sighat
 

        iatomstart = 1
        natomsp = natoms
        do ialp = iatomstart, iatomstart - 1 + natomsp
         rna(:) = ratom(:,ialp)
         indna = imass(ialp)

         do ineigh = 1, neigh_comn(ialp)
          mneigh = neigh_comm(ineigh,ialp)
 
          if (mneigh .ne. 0) then
           iatom = neigh_comj(1,ineigh,ialp)
           ibeta = neigh_comb(1,ineigh,ialp)
           r1(:) = ratom(:,iatom) + xl(:,ibeta)
           in1 = imass(iatom)

           jatom = neigh_comj(2,ineigh,ialp)
           jbeta = neigh_comb(2,ineigh,ialp)
           r2(:) = ratom(:,jatom) + xl(:,jbeta)
           in2 = imass(jatom)
           jneigh = neigh_back(iatom,mneigh)
           r21(:) = r2(:) - r1(:)
           y = sqrt(r21(1)*r21(1) + r21(2)*r21(2) + r21(3)*r21(3))
 
           r13(:) = r1(:) - rna(:)
           r23(:) = r2(:) - rna(:)
           distance_13 = sqrt(r13(1)**2 + r13(2)**2 + r13(3)**2)
           distance_23 = sqrt(r23(1)**2 + r23(2)**2 + r23(3)**2)

           if (y .lt. 1.0d-05) then
            sighat(1) = 0.0d0
            sighat(2) = 0.0d0
            sighat(3) = 1.0d0
            write (*,*) ' There is an error here in assemble_3c.f '
            write (*,*) ' r1 = r2!!!! BAD!!!! '
           else
            sighat(:) = r21(:)/y
           end if
 
           rnabc(:) = rna(:) - (r1(:) + r21(:)/2.0d0)
           x = sqrt(rnabc(1)**2 + rnabc(2)**2 + rnabc(3)**2)
 
           ! Find the unit vector in rnabc direction.
           if (x .lt. 1.0d-05) then
            rhat(1) = 0.0d0
            rhat(2) = 0.0d0
            rhat(3) = 0.0d0
           else
            rhat(:) = rnabc(:)/x
           end if
           cost = sighat(1)*rhat(1) + sighat(2)*rhat(2) + sighat(3)*rhat(3)
           call epsilon (rhat, sighat, eps)
!            if (igauss .eq. 0) then ! IF_DEF_ORDERN_END
             call trescentros (1, 0, isorpmax, in1, in2, indna, x, y, cost, eps, bcnax, nspecies)
           do inu = 1, num_orb(in2)
            do imu = 1, num_orb(in1)
             vna(imu,inu,mneigh,iatom) = vna(imu,inu,mneigh,iatom) + bcnax(imu,inu)*eq2
             vna(inu,imu,jneigh,jatom) = vna(imu,inu,mneigh,iatom)   !Symmetrize Hamiltonian (April 2018): jneigh is the  back_neigh:
            end do
           end do
! ****************************************************************************
! End loop over ialp and its common neighbors.
          end if
         end do ! do ineigh
        end do ! do ialp

! Format Statements
! ===========================================================================

        return
        end subroutine assemble_3c
