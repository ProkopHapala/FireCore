
! den2mesh.f90
! Program Description
! ===========================================================================
!       Project bands on the mesh.
!
!
!                 + X0 (iatom)
!                /|\     u1X = g1 - X0
! uX0 = X0- g0  / | \
!              /  |  + g1 (nearest grid point to iatom)
!             /   | /
!            /    |/
!           /     + Y0
!          +
!          g0 (origin of grid coords)
!
! g1 ... point where the origin of atomic mesh (sphere) is placed
!  vectors we need:
!    uX0 = X0 - g0
!    u1X = g1 - X0
!    r21 = Y0 - X0
!    u1Y = g1 - Y0 = g1 - Y0 - X0 + X0 = u1X - r21

subroutine mapAtomToGrid( u, X0, g1, index0 )
   use configuration
   use dimensions
   use interactions
   use neighbor_map
   use grid
   use density
   use charges
   use kpoints
   implicit none
! ============= Arguments
   real, dimension (3), intent(in) :: u
   real, dimension (3), intent(out):: X0
   real, dimension (3), intent(out):: g1
   integer, intent(out):: index0
! ============= Variables
   integer i0,j0,k0
   real, dimension (3) :: u1X
! ============= Body
   !u(:) = ratom(:,iatom) - g0(:)   ! vector between the atom and initial point of the grid
   !call mult3x1 (inva,u)            ! get n-vector
   ! round coefficients to get the position of the nearest grid point g1 to the iatom X1
   ! i,j,k can be positive or negative exceeding rmX (it means not centered in the unit cell)
   i0 = nint( u(1) )
   j0 = nint( u(2) )
   k0 = nint( u(3) )
   ! find the vector u1 between the iatom X1 and the nearest point g1
   u1X(1) = u(1) - real(i0)
   u1X(2) = u(2) - real(j0)
   u1X(3) = u(3) - real(k0)
   ! check if the nearest grid point is located within the unit cell of the grid coords
   ! if not, let's map it within
   !i0
   if (u(1) .lt. 0.0d0) then
      i0 = i0 + rm1*(int(abs(i0/rm1)) + 1)
   else
      i0 = i0 - rm1*int(i0/rm1)
   endif
   !j0
   if (u(2) .lt. 0.0d0) then
      j0 = j0 + rm2*(int(abs(j0/rm2)) + 1)
   else
      j0 = j0 - rm2*int(j0/rm2)
   endif
   !k0
   if (u(3) .lt. 0.0d0) then
      k0 = k0 + rm3*(int(abs(k0/rm3)) + 1)
   else
      k0 = k0 - rm3*int(k0/rm3)
   endif
   g1(:)  =             i0*elvec(1,:) +     j0*elvec(2,:) +     k0*elvec(3,:)   ! find the coordinates of the nearest point g1 witihin the grid coords
   X0(:)  = g1(:) + u1x(1)*elvec(1,:) + u1x(2)*elvec(2,:) + u1x(3)*elvec(3,:)   ! evaluate coordinates of the iatom in the grid coords
   index0 = 1 + (i0+emx1) + em1*(j0+emx2) + em1*em2*(k0+emx3)
end subroutine

!======================================================================
!======================================================================

subroutine getAtomBasis( in1, distX, dXr, psi1 )
    use configuration
    use dimensions
    use interactions
    use neighbor_map
    use grid
    use density
    use charges
    use kpoints
    implicit none
   ! ============= Arguments
   integer in1
   real distX
   real, dimension (3) :: dXr
   real, dimension(numorb_max), intent(out) :: psi1
   ! ============= Variables
   integer issh, l, lmu, imu
   real    psiR, dpsiR
   real, dimension (  5)  :: psiL
   real, dimension (3,5)  :: dpsiL
   ! ============= Body
   psi1 = 0.0d0
   imu = 1
   do issh = 1,nssh(in1)
      call getpsi(in1,issh,distX,psiR,dpsiR)   ! get radial part of wf.
      l = lssh(issh,in1)
      call getYlm(l,dXr,psiL,dpsiL)            ! get spherical harmonics part of wf.
      do lmu = 1, (2*l+1)
         psi1(imu) = psi1(imu) + psiL(lmu)*psiR
         imu = imu + 1
      enddo ! do lmu
   enddo ! do issh
end subroutine

!======================================================================
!======================================================================

subroutine project_dens( ewfaux, f_mul )
   use configuration
   use dimensions
   use interactions
   use neighbor_map
   use grid
   use density
   use charges
   use kpoints
   implicit none

! Argument Declaration and Description
! ===========================================================================
   !integer iband0,iband1
   real, dimension (nrm), intent (out) :: ewfaux
   ! real, target, dimension (:), allocatable :: ewfaux
   real, intent (in) :: f_mul

! Local Variable Declaration and Description
! ===========================================================================

   integer iatom, jatom, ineigh
   integer imu, inu
   integer in1, in2
   integer mbeta
   integer index, index0, ind
   !integer i, j, k
   integer imesh

   real distX, distY
   real dens

   real, dimension (3) :: r1, r2, r21
   real, dimension (3) :: u
   real, dimension (3) :: dXr, dYr
   real, dimension (3) :: X0
   real, dimension (3) :: g1
   real, dimension (3) :: u1X

!   real, dimension (nspec_max):: rcutoff_max
   real, dimension (numorb_max)   :: psi1,psi2
   !real, target, dimension (:), allocatable  :: rhotmp
   real, dimension (3,3)          :: amat
   real, dimension (3,3)          :: inva


! Procedure
! ===========================================================================
   amat = transpose(elvec)
   call inv3x3 (amat,inva)
   do iatom = 1, natoms
      in1   = imass(iatom)
      r1(:) = ratom(:,iatom)
      u(:) = ratom(:,iatom) - g0(:)   ! vector between the atom and initial point of the grid
      call mult3x1 (inva,u)            ! get n-vector
      call mapAtomToGrid( u, X0, g1, index0 )
      !write(*,*) "ia ", iatom, " i0 ", index0, " X0 ",  X0(:)," g1 ", g1(:)
      u1X(:) = g1(:) - X0(:)   ! vector pointing from g1 to X0
      ratom2g(:,iatom) = X0(:)  ! save iatom coord within the grid unit cell
      !index0 = 1 + (i0+emx1) + em1*(j0+emx2) + em1*em2*(k0+emx3)  ! find index of the gX point within the extended mesh
      ! Loop over the neighbors
      do ineigh = 1, neighn(iatom)
         jatom = neigh_j(ineigh,iatom)
         mbeta = neigh_b(ineigh,iatom)
         r2(:) = ratom(:,jatom) + xl(:,mbeta)
         in2 = imass(jatom)
         r21(:) = r2(:) - r1(:)
         ! Loop over points in the atomic mesh gP
         do imesh = 1, nam
            index  = index0 + am2rc(imesh)  ! restore index of the given mesh point gP within the extended mesh
            dXr(:) = ram2rc(:,imesh) + u1X(:)   ! evaluate the vector between the iatom and the mesh point gP
            dYr(:) = dXr(:) - r21(:)   ! evaluate the vector between the jatom and the mesh point gP
            distY  = sqrt(dYr(1)**2 + dYr(2)**2 + dYr(3)**2)  ! distance between the mesh point and jatom
            if (distY .lt. Rc_max) then
               distX = sqrt(dXr(1)**2 + dXr(2)**2 + dXr(3)**2)  ! distance between the mesh point and iatom
               call getAtomBasis( in1, distX, dXr, psi1 )
               call getAtomBasis( in2, distY, dYr, psi2 )
               ind = e2r(index) - 1  ! map the point from the extended mesh into the normal mesh
               dens = 0.0d0
               do inu = 1, num_orb(in1)
                  do imu = 1, num_orb(in2)
                     dens = dens + rho(inu,imu,ineigh,iatom)*psi1(inu)*psi2(imu)
                  enddo ! do inu
               enddo ! do imu
               ewfaux(ind) = ewfaux(ind) + dens * f_mul ! store variation of density at given point
            endif ! if (Rc_max)
         end do ! do imesh
      end do ! do ineigh
   end do ! do iatom
   !write (*,*) "project_orb vmin, vmax ", vmin, vmax
   return
end subroutine project_dens