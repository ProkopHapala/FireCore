
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
50          format (4i4, 2f12.6, 3f12.6)

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

