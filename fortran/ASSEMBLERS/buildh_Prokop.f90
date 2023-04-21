




subroutine buildh_prokop( ialp, icneigh, hmat, smat )
   ! ============ 2C

   ! ============ 3C
   do ialp = iatomstart, iatomstart - 1 + natomsp
      rna(:) = ratom(:,ialp)
      indna = imass(ialp)
      do icneigh = 1, neigh_comn(ialp)
         mneigh = neigh_comm(icneigh,ialp)
         if (mneigh .ne. 0) then
            call assemble_atom_neigh_3C( ialp, icneigh, hmat, smat )
         end if ! mneigh
      end do ! end do neigh_comn(ialp)
   end do  ! end do ialp
end subroutine assemble_atom_neigh()


! ===================================================
!                2C per atom-comb-neigh
!====================================================

subroutine assemble_atom_neigh_2C( iatom, ineigh, hmat, smat )
   use options
   implicit none
!======= arguments
   integer in1,in2

!======= Body

   matom = neigh_self(iatom)
   r1(:) = ratom(:,iatom)
   in1   = imass(iatom)

   mbeta = neigh_b(ineigh,iatom)
   jatom = neigh_j(ineigh,iatom)
   r2(:) = ratom(:,jatom) + xl(:,mbeta)
   in2   = imass(jatom)

   rcutoff_j = 0
   do imu = 1, nssh(in2)
      if (rcutoff(in2,imu) .gt. rcutoff_j) rcutoff_j = rcutoff(in2,imu)
   end do

   r21(:) = r2(:) - r1(:)
   y = sqrt(r21(1)*r21(1) + r21(2)*r21(2) + r21(3)*r21(3))

   if (y .lt. 1.0d-05) then
      sighat(1) = 0.0d0
      sighat(2) = 0.0d0
      sighat(3) = 1.0d0
   else
      sighat(:) = r21(:)/y
   end if

   call epsilon (r2, sighat, eps)
   call deps2cent (r1, r2, eps, deps)

   ! ========= FROM   assemble_olsxc_off.f90
   if (iatom .eq. jatom .and. mbeta .eq. 0) then
      ! Do nothing here - special case. Interaction already calculated in atm case.
   else
      call doscentros (6, 0, kforce, in1, in2, in2, y,     eps, deps, rhomx, rhompx)
      nnu = num_orb(in2)
      nmu = num_orb(in1)
      vxc(:nmu,:nnu,ineigh,iatom) =    vxc(:nmu,:nnu,ineigh,iatom) + rhomx(:nmu,:nnu,)
      do isorp = 1, nssh(in1)
         call doscentros (18, isorp, kforce, in1, in1, in2, y, eps, deps, rhomx, rhompx)
         dxn = (Qin(isorp,iatom) - Qneutral(isorp,in1))
         vxc_ca(:nmu,:nnu,ineigh,iatom) =  vxc_ca(:nmu,:nnu,ineigh,iatom) + rhomx(:nmu,:nnu,)*dxn
      end do
      ! the vxc_ontopr case, transfer of charge (McWeda)
      do isorp = 1, nssh(in2)
         call doscentros (19, isorp, kforce, in1, in2, in2, y,  eps, deps, rhomx, rhompx)
         dxn = (Qin(isorp,jatom) - Qneutral(isorp,in2))
         vxc_ca(:nmu,:nnu,ineigh,iatom) = vxc_ca(:nmu,:nnu,ineigh,iatom) + rhomx(:nmu,:nnu,)*dxn
      end do
      denmx(:nmu,:nnu) = rho_off  (:nmu,:nnu,ineigh,iatom)
      den1x(:nmu,:nnu) = rhoij_off(:nmu,:nnu,ineigh,iatom)
      sx   (:nmu,:nnu) = s_mat    (:nmu,:nnu,ineigh,iatom)
      call build_olsxc_off (in1, in2, den1x, denmx, sx, ineigh, iatom,  bcxcx)
   end if






   ! ========= FROM   assemble_ca_2c.f90
   call doscentros (9, 0, iforce, in1, in2, in2, y, eps, deps, dipx, dippx)
   nnu = num_orb(in2)
   nmu = num_orb(in1)
   dip(:nmu,:nnu,ineigh,iatom) = dipx(:nmu,:nnu)
   if (iforce .eq. 1) dipp(:,:nmu,:nnu,ineigh,iatom) = dippx(:,:nmu,:nnu)
   dq2 = 0.0d0
   do issh = 1, nssh(in2)
      dq2 = dq2 + (Qin(issh,jatom) - Qneutral(issh,in2))
   end do
   stn1 = 1.0d0
   stn2 = 0.0d0
   emnpl = 0.0d0
   if (y .gt. 1.0d-4) then
      icount_sav = 0
      do issh = 1, nssh(in1)
         jcount_sav = 0
         do jssh = 1, nssh(in1)
            rend = rcutoff_i + rcutoff_j
            call smoother (y, rend, smt_elect, stn_temp1, dstn_temp)
            stn_temp2 = 1.0d0 - stn_temp1
            do inu = 1, lssh(issh,in1)*2 + 1
               icount = icount_sav + inu
               do imu = 1, lssh(jssh,in1)*2 + 1
                  jcount = jcount_sav + imu
                  stn1(icount,jcount) = stn_temp1
                  stn2(icount,jcount) = stn_temp2
               end do
            end do
            jcount_sav = jcount
         end do
         icount_sav = icount
      end do
      nnu = num_orb(in1)
      nmu = num_orb(in1)
      emnpl(:nmu,:nnu,) = (s_mat(:nmu,:nnu,,matom,iatom)/y)*dq2
      ewaldsr(:nmu,:nnu,,matom,iatom) =   ewaldsr(:nmu,:nnu,,matom,iatom) + emnpl(:nmu,:nnu,)*eq2
   end if  ! End if y .gt. 1.0d-4
   nnu = num_orb(in3)
   nmu = num_orb(in1)
   kforce = 0
   interaction = 4
   in3 = in1
   do isorp = 1, nssh(in2)
      call doscentros (4, isorp, 0, in1, in2, in1, y,  eps, deps, bccax, bccapx)
      dxn = (Qin(isorp,jatom) - Qneutral(isorp,in2))
      bcca(:nmu,:nnu) = bcca(:nmu,:nnu) + bccax(:nmu,:nnu)*dxn
   end do
   vca(:nmu,:nnu,matom,iatom) = vca(:nmu,:nnu,matom,iatom)  + (stn1(:nmu,:nnu)*bcca(:nmu,:nnu) + stn2(:nmu,:nnu)*emnpl(:nmu,:nnu))*eq2
   if (y .gt. 1.0d-4) then
      nnu = num_orb(in2)
      nmu = num_orb(in1)
      ewaldsr(:nmu,:nnu,ineigh,iatom) = ewaldsr(:nmu,:nnu,ineigh,iatom)  + ( ( (s_mat(:nmu,:nnu,ineigh,iatom)/(2.0d0*y))    +  (dip(:nmu,:nnu,ineigh,iatom)/(y*y))  )*dq1 + (  (dq2*s_mat(:nmu,:nnu,ineigh,iatom)/(2.0d0*y)) - (dip(:nmu,:nnu,ineigh,iatom)/(y*y))   )*dq2  )*eq2
   end if
   if (iatom .eq. jatom .and. mbeta .eq. 0) then
      ! Do nothing here - special case. Interaction already calculated in atm case.
   else
      bcca = 0.0d0
      nnu = num_orb(in3)
      nmu = num_orb(in1)
      do isorp = 1, nssh(in1)
         call doscentros (2, isorp, kforce, in1, in1, in2, y, eps, deps, bccax, bccapx)
         dxn = (Qin(isorp,iatom) - Qneutral(isorp,in1))
         bcca(:nmu,:nnu) = bcca(:nmu,:nnu) + dxn*bccax(:nmu,:nnu)
      end do
      do isorp = 1, nssh(in2)
         call doscentros (3, isorp, kforce, in1, in2, in2, y,  eps, deps, bccax, bccapx)
         dxn = (Qin(isorp,jatom) - Qneutral(isorp,in2))
         bcca(:nmu,:nnu) = bcca(:nmu,:nnu) + dxn*bccax(:nmu,:nnu)
      end do
      vca(:nmu,:nnu,ineigh,iatom) =   vca(:nmu,:nnu,ineigh,iatom) + bcca(:nmu,:nnu)*eq2
   end if
   
end subroutine assemble_atom_neigh_2C()


! ===================================================
!                3C per atom-comb-neigh
!====================================================

subroutine assemble_atom_neigh_3C( ialp, icneigh, hmat, smat )
   use options
   implicit none
!======= arguments
   integer in1,in2

!======= Body

   if (Kscf .eq. 1) then
      call cl_value (indna, cl)
   end if


   rna(:) = ratom(:,ialp)
   indna  = imass(ialp)
   do ineigh = 1, neigh_comn(ialp)
      mneigh = neigh_comm(ineigh,ialp)
      if (mneigh .ne. 0) then
         iatom = neigh_comj(1,icneigh,ialp)
         ibeta = neigh_comb(1,icneigh,ialp)
         r1(:) = ratom(:,iatom) + xl(:,ibeta)
         in1 = imass(iatom)

         jatom = neigh_comj(2,icneigh,ialp)
         jbeta = neigh_comb(2,icneigh,ialp)
         r2(:) = ratom(:,jatom) + xl(:,jbeta)
         in2 = imass(jatom)
         jneigh = neigh_back(iatom,mneigh)

         ! Find r21 = vector pointing from r1 to r2, the two ends of the bondcharge. This gives us the distance dbc (or y value in the 2D grid).
         r21(:) = r2(:) - r1(:)
         y = sqrt(r21(1)*r21(1) + r21(2)*r21(2) + r21(3)*r21(3))

         ! Find the unit vector in sigma direction.
         if (y .lt. 1.0d-05) then
            sighat(1) = 0.0d0
            sighat(2) = 0.0d0
            sighat(3) = 1.0d0
            write (*,*) ' There is an error here in assemble_3c.f '
            write (*,*) ' r1 = r2!!!! BAD!!!! '
         else
            sighat(:) = r21(:)/y
         end if

         ! Find rnabc = vector pointing from center of bondcharge to the neutral atom. ! The center of the bondcharge is at r1 + r21/2.  This gives us the distance dnabc (x value in 2D grid).
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


         do isorp = 1, nssh(indna)
            ! ========== FROM average_rho_ca.f90
            call trescentros ( 3, isorp, isorpmax, in1, in2, indna, x, y, cost, eps, rhomx, nspecies)
            nnu = num_orb(in2)
            nmu = num_orb(in1)
            rho_off(:nmu,:nnu,mneigh,iatom) = rho_off(:nmu,:nnu,mneigh,iatom) + rhomx(:nmu,:nnu,)*Qin(isorp,ialp)
            rho_off(:nmu,:nnu,jneigh,jatom) = rho_off(:nmu,:nnu,mneigh,iatom)  !Symmetrize:
            call trescentrosS(             isorp, isorpmax, in1, in2, indna, x, y, cost, eps, rhomm, nspecies)
            rhom_3c(:nmu,:nnu,mneigh,iatom) = rhom_3c(:nmu,:nnu,mneigh,iatom) + rhomm(:nmu,:nnu,)*Qin(isorp,ialp)
            rhom_3c(:nmu,:nnu,jneigh,jatom) = rhom_3c(:nmu,:nnu,mneigh,iatom) !Symmetrize:
            ! ========== END average_rho_ca.f90
            if (Kscf .eq. 1) then
               ! ========= FROM assemble_3c.f90
               call trescentros (1, 0, isorpmax, in1, in2, indna, x, y, cost, eps, bcnax, nspecies)
               vna(:nmu,:nnu,mneigh,iatom) = vna(:nmu,:nnu,mneigh,iatom) + bcnax(:nmu,:nnu)*eq2
               ! ========= END  assemble_3c.f90
               ! ========= FROM assemble_3c_PP.f90
               r31(:) = rna(:) - r1(:)
               r32(:) = rna(:) - r2(:)
               m31 = mpairnay(iatom, ialp, r31)
               m32 = mpairnay(jatom, ialp, r32)
               bcnlx(imu,inu) = 0.0d0
               do ncc = 1, num_orbPP(indna)
                  bcnlx(:nmu,:nnu) = bcnlx(:nmu,:nnu) + cl(ncc)*sVNL(:nmu,ncc,m31,iatom)*sVNL(:nnu,ncc,m32,jatom)
               end do
               vnl(:nmu,:nnu,mneigh,iatom) = vnl(:nmu,:nnu,mneigh,iatom) + bcnlx(:nmu,:nnu)
               ! ========= END assemble_3c_PP.f90
            end if ! if (Kscf .eq. 1)
         end do ! do isorp

         ! =========================================
         ! ========= FROM   assemble_ca_3c.f90
         distance_13 = sqrt((rna(1) - r1(1))**2 + (rna(2) - r1(2))**2 + (rna(3) - r1(3))**2)
         distance_23 = sqrt((rna(1) - r2(1))**2 + (rna(2) - r2(2))**2 + (rna(3) - r2(3))**2)
         icount_sav = 0
         do issh = 1, nssh(in1)
            jcount_sav = 0
            do jssh = 1, nssh(in2)
               rend1 = rcutoff_i + rcutoff_ialp
               rend2 = rcutoff_j + rcutoff_ialp
               call smoother (distance_13, rend1, smt_elect, stn_temp1, dstn_temp)
               call smoother (distance_23, rend2, smt_elect, stn_temp2, dstn_temp)
               stn_temp1 = stn_temp1*stn_temp2
               stn_temp2 = 1.0d0 - stn_temp1
               do inu = 1, lssh(issh,in1)*2 + 1
                  icount = icount_sav + inu
                  do imu = 1, lssh(jssh,in2)*2 + 1
                     jcount = jcount_sav + imu
                     stn1(icount,jcount) = stn_temp1
                     stn2(icount,jcount) = stn_temp2
                  end do
               end do
               jcount_sav = jcount
            end do ! jssh
            icount_sav = icount
         end do ! issh
         do inu = 1, num_orb(in2)
            do imu = 1, num_orb(in1)
               sterm = dq3*s_mat(imu,inu,mneigh,iatom)/2.0d0
               dterm = dq3*dip  (imu,inu,mneigh,iatom)/y
               emnpl  (imu,inu)              = (sterm - dterm)/distance_13 + (sterm + dterm)/distance_23
               ewaldsr(imu,inu,mneigh,iatom) = ewaldsr(imu,inu,mneigh,iatom) + emnpl(imu,inu)*eq2
               ewaldsr(inu,imu,jneigh,jatom) = ewaldsr(imu,inu,mneigh,iatom)
            end do
         end do
         nnu = num_orb(in2)
         nmu = num_orb(in1)
         bcca = 0.0d0
         do isorp = 1, nssh(indna)
            interaction = 1
            call trescentros (interaction, isorp, isorpmax, in1, in2, indna, x, y, cost, eps, bccax, nspecies)
            dxn = (Qin(isorp,ialp) - Qneutral(isorp,indna))
            bcca(:nmu,:nnu) = bcca(:nmu,:nnu) + bccax(:nmu,:nnu)*dxn
         end do
         vca(:nmu,:nnu,mneigh,iatom) = vca(:nmu,:nnu,mneigh,iatom) + ( stn1(:nmu,:nnu)*bcca(:nmu,:nnu) + stn2(:nmu,:nnu)*emnpl(:nmu,:nnu) )*eq2
         vca(:nmu,:nnu,jneigh,jatom) = vca(:nmu,:nnu,mneigh,iatom)
         ! ========= END   assemble_ca_3c.f90


      end if
   end do ! end do neigh_comn(ialp)
end do  ! end do ialp

end subroutine assemble_atom_neigh()





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

   real distance

   real, dimension (numorb_max, numorb_max) :: htemp
   real, dimension (numorb_max, numorb_max) :: stemp

   real, dimension (3) :: dvec

   integer issh
   integer numorb
   integer jatom0
   integer ineigh0
   integer mbeta0

! ===========================================================================
! Procedure
! ===========================================================================


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

   if (Kscf .eq. 1) then
      call reallocate_neigh ! (nprocs, my_proc, iordern, itheory, itheory_xc, igauss, icluster, ivdw, iwrthampiece, iwrtatom, igrid)
      call neighbors   ! (nprocs, my_proc, iordern, icluster, iwrtneigh, ivdw)
      call neighborsPP ! (nprocs, my_proc, iordern, icluster,iwrtneigh)
      call num_neigh_tot ! (numorb_max)
      call backnay ()
      call neighbors_pairs ! (icluster)
      call common_neighbors   ! (nprocs, my_proc, iordern, iwrtneigh_com)
      call common_neighborsPP ! (nprocs, my_proc, iordern,  iwrtneigh_com, icluster)
      call check_neighbors()
   end if ! end if (Kscf .eq. 1)

   ! ========================== BIG LOOP


   iatomstart = 1
   natomsp = natoms

   do iatom = iatomstart, iatomstart - 1 + natomsp
      in1 = imass(iatom)




      ! 1) =========== From   assemble_olsxc_1c()
      call unocentros (in1, iatom, iforce, itheory, ixc, exc_1c, muexc_1c,  dccexc_1c, mu1xc)
      etotxc_1c = etotxc_1c + dccexc_1c
      do imu = 1, num_orb(in1)
         do inu = 1, num_orb(in1)
            vxc_1c(imu,inu,matom,iatom) =  vxc_1c(imu,inu,matom,iatom) + mu1xc(imu,inu)
         end do
      end do


      do ineigh = 1, neighn(iatom)
         mbeta = neigh_b(ineigh,iatom)
         jatom = neigh_j(ineigh,iatom)
         in2 = imass(jatom)
         r2(:) = ratom(:,jatom) + xl(:,mbeta)

         ! Find the unit vector in sigma direction.
         if (y .lt. 1.0d-05) then
            sighat(1) = 0.0d0
            sighat(2) = 0.0d0
            sighat(3) = 1.0d0
         else
            sighat(:) = r21(:)/y
         end if
         call epsilon  (r2,     sighat,  eps)
         call deps2cent(r1, r2, eps,    deps)
         r21(:) = r2(:) - r1(:)
         y = sqrt(r21(1)*r21(1) + r21(2)*r21(2) + r21(3)*r21(3))


         ! 2) =========== From assemble_sVNL()   ! Call doscentros and get 2-center "overlap" interactions of VNL.  ! That is we get <phi_i | VNL_j>. We treat it just like an overlap.
         isorp = 0
         interaction = 5
         call doscentrosPP (interaction, isorp, y, eps, deps, iforce, in1, in2, sVNLx, spVNLx)
         do inu = 1, num_orbPP(in2)
            do imu = 1, num_orb(in1)
               if (ineigh .ne. matom) then
                  sVNL(imu,inu,ineigh,iatom) = sVNLx(imu,inu)
                  if (iforce .eq. 1) then
                     spVNL(:,imu,inu,ineigh,iatom) = spVNLx(:,imu,inu)
                  end if
               else
                  sVNL(imu,inu,ineigh,iatom) = sVNLx(imu,inu)
                  spVNL(:,imu,inu,ineigh,iatom) = 0.0d0
               end if
            end do
         end do



! =======================  FROM  subroutine assemble_2c()
! ================== Kinetic Energy & Overlap
         isorp = 0
         interaction = 1   ! Overlap
         in3 = in2
         call doscentros (interaction, isorp, iforce, in1, in2, in3, y,  eps, deps, sx, spx)
         !call debug_writeIntegralSet( "intS.log", interaction, isorp, in1, in2 )
         isorp = 0
         interaction = 13   ! Kinetic
         in3 = in2
         call doscentros (interaction, isorp, iforce, in1, in2, in3, y, eps, deps, tx, tpx)
! Write s and t to appropriate arrays
         do inu = 1, num_orb(in2)
            do imu = 1, num_orb(in1)
               s_mat(imu,inu,ineigh,iatom) = sx(imu,inu)
               t_mat(imu,inu,ineigh,iatom) = tx(imu,inu)
               if (iforce .eq. 1) then
                  sp_mat(:,imu,inu,ineigh,iatom) = spx(:,imu,inu)
                  tp_mat(:,imu,inu,ineigh,iatom) = tpx(:,imu,inu)
                  !write(*,*) "i,j,inu,imu,spx,tpx ", iatom,jatom,inu,imu,  spx(:,imu,inu), tpx(:,imu,inu)
               end if
            end do
         end do
! ================== CALL DOSCENTROS AND GET VNA FOR ATOM CASE ! The vna 2 centers are: ontop (L), ontop (R), and atm. ! First, do vna_atom case. Here we compute <i | v(j) | i> matrix elements.
         isorp = 0
         kforce = 1                             ! don't calculate forces here
         interaction = 4
         in3 = in1
         call doscentros (interaction, isorp, kforce, in1, in2, in3, y,  eps, deps, bcnax, bcnapx)
         do inu = 1, num_orb(in3)
            do imu = 1, num_orb(in1)
               vna(imu,inu,matom,iatom) = vna(imu,inu,matom,iatom) + bcnax(imu,inu)*eq2
            end do
         end do
! ================== VNA FOR ONTOP CASE
         if (iatom .eq. jatom .and. mbeta .eq. 0) then
         else  ! Do nothing here - special case. Interaction already calculated in atm case.
            isorp = 0
            interaction = 2
            in3 = in2
            call doscentros (interaction, isorp, kforce, in1, in1, in3,  y, eps, deps, bcnax, bcnapx)
            do inu = 1, num_orb(in3)
               do imu = 1, num_orb(in1)
                  bcna(imu,inu) = bcnax(imu,inu)
               end do
            end do
! For the vna_ontopr case, the potential is in the second atom (jatom): ! Neutral atom piece
            isorp = 0
            interaction = 3
            in3 = in2
            call doscentros (interaction, isorp, kforce, in1, in2, in3,  y, eps, deps, bcnax, bcnapx)
            do inu = 1, num_orb(in3)
               do imu = 1, num_orb(in1)
                  bcna(imu,inu) = bcna(imu,inu) + bcnax(imu,inu)
                  vna(imu,inu,ineigh,iatom) = vna(imu,inu,ineigh,iatom) + bcna(imu,inu)*eq2
               end do
            end do

         end if ! End if for r1 .ne. r2 case


         !!!!! MISSING ToDo !!!!! MISSING ToDo !!!!! MISSING ToDo
         if (itheory .eq. 1) then
            call average_ca_rho()
         else
            call average_rho()
         endif
         !!!!! MISSING ToDo !!!!! MISSING ToDo !!!!! MISSING ToDo







         kneigh = nPPx_map(ineigh,iatom)

! Here we add in all of the charge interactions for Kohn-Sham method .
         if (itheory .eq. 3) then
            do inu = 1, num_orb(in2)
               do imu = 1, num_orb(in1)



                  ! =============== From assemble_2c_PP()
                  ! ASSEMBLE VNL ATM CASE  <phi_i|Psi_j><Psi_j|phi_i>
                  call cl_value (in2, cl)
                  do ncc = 1, num_orbPP(in2)
                     PPx(imu,inu) = PPx(imu,inu) + cl(ncc)*sVNL(imu,ncc,ineigh,iatom)*sVNL(inu,ncc,ineigh,iatom)
                  end do
                  vnl(imu,inu,matom,iatom) = vnl(imu,inu,matom,iatom) + PPx(imu,inu)
                  ! ASSEMBLE VNL ONTOP LEFT CASE   <phi_i|Psi_i><Psi_i|phi_j>
                  if (iatom .eq. jatom .and. mbeta .eq. 0) then      ! Case 1. PP is iatom.  <i | VNL(i) |j>.
                     call cl_value (in1, cl)
                     jneigh = nPPx_point(ineigh,iatom)
                     mneigh_self = nPP_self(iatom)
                     PPx(imu,inu) = 0.0d0
                     do ncc = 1, num_orbPP(in1)
                        PPx(imu,inu) = PPx(imu,inu) + cl(ncc)*sVNL(imu,ncc,mneigh_self,iatom) *sVNL(inu,ncc,jneigh,jatom)
                     end do ! do ncc
                     vnl(imu,inu,kneigh,iatom) = vnl(imu,inu,kneigh,iatom) + PPx(imu,inu)
                     ! ASSEMBLE VNL ONTOP RIGHT CASE   <phi_i|Psi_j><Psi_j|phi_j>
                     if (iatom .eq. jatom .and. mbeta .eq. 0) then
                     else ! if(iatom .eq. jatom)
                        call cl_value (in2, cl)
                        mneigh_self = nPP_self(jatom)
                        PPx(imu,inu) = 0.0d0
                        do ncc = 1, num_orbPP(in2)
                           PPx(imu,inu) = PPx(imu,inu) + cl(ncc)*sVNL(imu,ncc,ineigh,iatom) *sVNL(inu,ncc,mneigh_self,jatom)
                        end do
                        vnl(imu,inu,kneigh,iatom) = vnl(imu,inu,kneigh,iatom) + PPx(imu,inu)



                        !  ToDo : Later we will get-rid off intermediate arrays   t_mat[], vna[], vxc[], vca[]
                        h_mat(imu,inu,ineigh,iatom) =  h_mat(imu,inu,ineigh,iatom) + vca(imu,inu,ineigh,iatom) + vxc_ca(imu,inu,ineigh,iatom) + ewaldlr(imu,inu,ineigh,iatom) - ewaldsr(imu,inu,ineigh,iatom)      !  if (itheory .eq. 1 .or. itheory .eq. 2) then
                        h_mat(imu,inu,ineigh,iatom) =  t_mat(imu,inu,ineigh,iatom) + vna(imu,inu,ineigh,iatom) + vxc(imu,inu,ineigh,iatom)    + vca(imu,inu,ineigh,iatom)      ! (itheory .eq. 3) then
                        h_mat(imu,inu,ineigh,iatom) =  t_mat(imu,inu,ineigh,iatom) + vna(imu,inu,ineigh,iatom) + vxc(imu,inu,ineigh,iatom)    + vxc_1c(imu,inu,ineigh,iatom)   ! else


                     end do
                  end do




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
                              do inu = 1, num_orb(in2)
                                 do imu = 1, num_orb(in1)
                                    h_mat(imu,inu,ineigh,iatom) = 0.0d0
                                    t_mat(imu,inu,ineigh,iatom) = 0.0d0
                                    s_mat(imu,inu,ineigh,iatom) = 0.0d0
                                    vna(imu,inu,ineigh,iatom) = 0.0d0
                                    vxc(imu,inu,ineigh,iatom) = 0.0d0
                                    if (itheory .eq. 1 .or. itheory .eq. 2) then
                                       ewaldlr(imu,inu,ineigh,iatom) = 0.0d0
                                       ewaldsr(imu,inu,ineigh,iatom) = 0.0d0
!               ewaldqmmm(imu,inu,ineigh,iatom) = 0.0d0    ! IF_DEF_QMMM_END
                                    end if
                                 end do
                              end do
                           end if
                        end if
                     end do ! do ineigh
                  end do ! do iatom

! Format Statements
! ===========================================================================
50                format (4i4, 2f12.6, 3f12.6)

                  return
               end subroutine buildh
























               FROM    assemble_mcweda ()
               if ((itheory .eq. 1 .or. itheory .eq. 2) .and. Kscf .eq. 1) then
                  call get_ewald ! (nprocs, my_proc, kforce, icluster, itheory, iordern)
               end if
               if(itheory_xc .eq. 2 ) then
                  call assemble_olsxc_1c ! (natoms, itheory, iforce)
               endif
               call timer_start_i(8)
               if (Kscf .eq. 1) then
                  call assemble_sVNL()
                  call assemble_2c()
                  call assemble_2c_PP()
                  call timer_stop_i(20)
               end if
               if (itheory_xc .eq. 1 ) then
                  !write (*,*) ' Assemble SN-xc exchange-correlation interactions. '
                  if (itheory .eq. 1) then
                     call average_ca_rho()
                  else
                     call average_rho()
                  endif
                  call assemble_snxc_on ( uxcdcc_sn )
                  call timer_start_i(23)
                  call assemble_snxc_off()
                  call timer_stop_i(23)
               end if ! if (itheory_xc = 1)
               if (itheory_xc .eq. 2 ) then
                  if (itheory .eq. 1) then
                     call average_ca_rho ! (nprocs, Kscf, iforce, iordern, igauss)
                  else
                     call average_rho ! (nprocs, Kscf, iforce, iordern, igauss)
                  endif
                  call assemble_olsxc_on (uxcdcc_ols) ! (natoms, nprocs, my_proc, iordern, itheory, uxcdcc_ols)
                  call assemble_olsxc_off ! (nprocs, my_proc, iordern, itheory)
               end if ! if (itheory_xc = 2)
               if (itheory .eq. 1) then
                  call assemble_ca_2c()
               endif
               if (Kscf .eq. 1) then
                  call assemble_3c    ! (nprocs, iordern, igauss, itheory_xc)
                  call assemble_3c_PP ! (nprocs, iordern)
               end if
               if (itheory .eq. 1) then
                  call assemble_ca_3c()
                  call assemble_lr   ()

