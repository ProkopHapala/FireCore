

! ===================================================
! ===================================================

subroutine write_grid_info() !(icluster)
    use configuration
    use grid
    implicit none
! ======= Variables
! ======= Body
    write(*,*) "!================ write_grid_info() "
    write(*,*) "regular  mesh size(:) ",   rm1,rm2,rm3, " tot ", nrm
    write(*,*) "extenden mesh size(:) ",   em1,em2,em3, " tot ", nem
    write(*,*) "atomic   mesh size(:) ",   am1,am2,am3, " tot ", nam
    write(*,*) "FDM      mesh size(:) ",   mfd1,mfd2,mfd3, " tot ", nmfd
    write(*,*) "dvol ", dvol, "dvol ", volume, " drmax ", drmax
    write(*,*) "g0 ", g0(:)
    write(*,*) "   avec(:,:) "
    write(*,*) a1vec(:)
    write(*,*) a2vec(:)
    write(*,*) a3vec(:)
    write(*,*) "   rlvec(:,:) "
    write(*,*) rlvec(1,:)
    write(*,*) rlvec(2,:)
    write(*,*) rlvec(3,:)
    write(*,*) "   elvec(:,:)    "
    write(*,*) elvec(1,:)
    write(*,*) elvec(2,:)
    write(*,*) elvec(3,:)
end subroutine

! ===================================================
! ===================================================

subroutine invertCell( a1vec,a2vec,a3vec, rlvec ) !(icluster)
    implicit none
    ! ======= Parameters
    real, dimension (3), intent(in) :: a1vec
    real, dimension (3), intent(in) :: a2vec
    real, dimension (3), intent(in) :: a3vec
    real, dimension (3,3), intent(out) :: rlvec
    ! ======= Variables
    integer i
    real denom
    real, dimension (3)     :: acrossb
    real, dimension (3)     :: bcrossc
    real, dimension (3)     :: ccrossa
    ! ======= Body
    call crossx(a2vec,a3vec,bcrossc)
    call crossx(a3vec,a1vec,ccrossa)
    call crossx(a1vec,a2vec,acrossb)
    denom=a1vec(1)*bcrossc(1)+a1vec(2)*bcrossc(2)+a1vec(3)*bcrossc(3)
    do i=1,3
        rlvec(1,i)=bcrossc(i)/denom
        rlvec(2,i)=ccrossa(i)/denom
        rlvec(3,i)=acrossb(i)/denom
    enddo
end subroutine

! ===================================================
! ===================================================

subroutine gridSizeFromRcut() !(icluster)
    use options, only: icluster
    use grid
    use constants_fireball
    use interactions
    use configuration
    use dimensions
    implicit none
    ! ======= Variables
    integer i
    real, dimension (3)     :: cvec
    ! ======= Body
    do i = 1,3
        cvec(i) = sqrt (rlvec(i,1)**2 + rlvec(i,2)**2 + rlvec(i,3)**2)        ! calc the size of the unit cell along the i-axis
        ngrid(i)   = int (2 * sqrt(Ecut) / (cvec(i)*abohr) + 1)                    ! estimate the number of points along the i-axis
        do     ! iterate until the number of points is multiplier of 2
        if ( mod( ngrid(i), 2 ) .eq. 0 ) exit     ! number of points must be multiply of 2!!
        ngrid(i) = ngrid(i) + 1
        enddo
        ervec(i,:) = (rlvec(i,:) * ngrid(i)) / (2.0d0*pi)
        cvec(i) = cvec(i) / ngrid(i)    ! the elementary distance at given axis Should we multiply by np  here as well? (In order to compute drmax a few lines later)
    enddo ! do i
    ! save max dr
    drmax = 0.0d0
    do i = 1,3
        if (drmax .gt. cvec(i)) drmax = cvec(i)
    enddo
end subroutine

! ===================================================
! ===================================================

subroutine sizeGrid !(icluster)
    use options, only: icluster
    use grid
    use constants_fireball
    use interactions
    use configuration
    use dimensions
    implicit none
! ======= Variables
    integer i
    real, dimension (3)     :: cvec
! ======= Body

 ! calculate volume of elementar unit cell
cvec(1)= a1vec(2)*a2vec(3) - a1vec(3)*a2vec(2)
cvec(2)= a1vec(3)*a2vec(1) - a1vec(1)*a2vec(3)
cvec(3)= a1vec(1)*a2vec(2) - a1vec(2)*a2vec(1)
volume = abs (a3vec(1)*cvec(1) + a3vec(2)*cvec(2) + a3vec(3)*cvec(3))
write(*,*) "volume ", volume

call invertCell( a1vec,a2vec,a3vec, rlvec )
rlvec(:,:) = rlvec(:,:)*(2.0d0*pi)
call  gridSizeFromRcut()
call invertCell( ervec(1,:),ervec(2,:),ervec(3,:), elvec )

write(*,*) "   elvec(:,:)    "
write(*,*) elvec(1,:)
write(*,*) elvec(2,:)
write(*,*) elvec(3,:)



!call crossx(a2vec,a3vec,bcrossc)
!call crossx(a3vec,a1vec,ccrossa)
!call crossx(a1vec,a2vec,acrossb)
!denom=a1vec(1)*bcrossc(1)+a1vec(2)*bcrossc(2)+a1vec(3)*bcrossc(3)
!do i=1,3
!   rlvec(1,i)=2.0d0*pi*bcrossc(i)/denom
!   rlvec(2,i)=2.0d0*pi*ccrossa(i)/denom
!   rlvec(3,i)=2.0d0*pi*acrossb(i)/denom
!enddo

!er1vec(:) = ervec(1,:)
!er2vec(:) = ervec(2,:)
!er3vec(:) = ervec(3,:)
! Find direct elementary vectors
!call crossx(er2vec,er3vec,bcrossc)
!call crossx(er3vec,er1vec,ccrossa)
!call crossx(er1vec,er2vec,acrossb)
!denom= er1vec(1)*bcrossc(1) + er1vec(2)*bcrossc(2) + er1vec(3)*bcrossc(3)
!do i=1,3
!   elvec(1,i)=bcrossc(i)/denom
!   elvec(2,i)=ccrossa(i)/denom
!   elvec(3,i)=acrossb(i)/denom
!enddo

rm1 = ngrid(1)  ! store the division of the regular mesh along the axis
rm2 = ngrid(2)
rm3 = ngrid(3)
nrm = rm1 * rm2 * rm3   ! total number of points on the regular mesh
! elementary volume
dvol = abs(elvec(1,1)*(elvec(2,2)*elvec(3,3)-elvec(2,3)*elvec(3,2))    &
&           +elvec(1,2)*(elvec(2,3)*elvec(3,1)-elvec(2,1)*elvec(3,3))    &
&           +elvec(1,3)*(elvec(2,1)*elvec(3,2)-elvec(2,2)*elvec(3,1)))
end subroutine

! ===================================================
! ===================================================

subroutine center_cell !(icluster)
    use options, only: icluster
    use grid
    use constants_fireball
    use interactions
    use configuration
    use dimensions
    implicit none
! ======= Variables
   real, parameter :: droff = 0.2d0
   integer i,j,ix,iatom, ispec
   real xmax,xmin
   !   real, parameter :: droff = 1.5d0
    real, dimension (3,3)     :: avec
    real, dimension (3)       :: cmass
! ======= Body

    ! ---- find Rc_max
    Rc_max = 0.0d0
    do ispec = 1, nspecies
    do i = 1, nssh(ispec)
        if (rcutoff(ispec,i) .gt. Rc_max) Rc_max = rcutoff(ispec,i)
    end do
    end do
    write(*,*) "Rc_max ", Rc_max

   if (icluster .eq. 1) then   !   non-periodic boundary conditions
    avec = 0.0d0
    do ix = 1,3   ! find max and min position of atom in ix-axis direction
       xmax = -8000.0d0
       xmin = 8000.0d0
       do iatom = 1,natoms
          if (xmax .lt. ratom(ix,iatom)) xmax = ratom(ix,iatom)
          if (xmin .gt. ratom(ix,iatom)) xmin = ratom(ix,iatom)
       enddo ! do iatom
       avec(ix,ix) = xmax - xmin + 2*Rc_max + 2*droff     ! define lattice vector
       if (ifixg0 .eq. 0) g0(ix) = xmin - Rc_max - droff  ! define intial point
    enddo ! do ix
    a1vec(:) = 0.0d0        ! setup lattice vector
    a1vec(1) = avec(1,1)
    a2vec(:) = 0.0d0
    a2vec(2) = avec(2,2)
    a3vec(:) = 0.0d0
    a3vec(3) = avec(3,3)
    do i = 1,natoms    ! find center of mass of atoms in the unit cell
       do j = 1,3
          cmass(j) = cmass(j) + ratom(j,i)
       enddo ! do j
    enddo ! do i
    cmass(:) = cmass(:)/float(natoms)
 else
    avec(1,:) = a1vec(:)
    avec(2,:) = a2vec(:)
    avec(3,:) = a3vec(:)
    do i = 1,natoms    ! find center of mass of atoms in the unit cell
       do j = 1,3
          cmass(j) = cmass(j) + ratom(j,i)
       enddo ! do j
    enddo ! do i
    cmass(:) = cmass(:)/float(natoms)
    if (ifixg0 .eq. 0) then
       do i = 1,3      ! find  initial point of the mesh. defined through the center of mass
          g0(i) = cmass(i) - 0.5d0*(a1vec(i) + a2vec(i) + a3vec(i))
       enddo
    endif ! if (ifixg0 .eq. 1)
 endif ! if (icluster .eq. 1)
 write (*,*) ' ---------------------------------'
 write (*,*) '          Lattice vector          '
 write (*,*) ' ---------------------------------'
 do i = 1,3
    write (*,*) (avec(i,j),j=1,3)
 enddo
 write (*,*)  'atomic units'
 do i = 1,3
    write (*,*) (avec(i,j)/abohr,j=1,3)
 enddo
end subroutine

! ===================================================
! ===================================================

subroutine initgrid_new !(icluster)
   use options, only: icluster
   use grid
   use constants_fireball
   use interactions
   use configuration
   use dimensions
   implicit none

! Argument Declaration and Description
! ===========================================================================
! Input
!   integer, intent (in)                     :: icluster

! Output

! Local Parameters and Data Declaration
! ===========================================================================

! Local Variable Declaration and Description
! ===========================================================================

   integer n1
   integer n2
   integer n3
   integer ne1
   integer ne2
   integer i
   integer j
   integer k
   integer ix
   integer index
   integer index0
   integer ispec
   integer ipoint
   integer iatom
   integer ii
   integer jj
   integer kk

   real r1
   real r2
   real r3
   real dist
   real xmin
   real xmax
   real denom
   !real drmax
   real coxm

   real, dimension (3)       :: cmass
   real, dimension (3)       :: cvec
   real, dimension (3)       :: xvec
   real, dimension (3,3)     :: avec
   real, dimension (3,3)     :: rvec
   integer, dimension (3)    :: ipiv
   integer, dimension (3)    :: np
   !real, dimension (3)      :: er1vec
   !real, dimension (3)      :: er2vec
   !real, dimension (3)      :: er3vec
   !real, dimension (3)      :: acrossb
   !real, dimension (3)      :: bcrossc
   !real, dimension (3)      :: ccrossa


! Allocate arrays
! ===========================================================================

! Procedure
! ===========================================================================
!   write (*,*) "DEBUG initgrid "

    call center_cell() 
    call sizeGrid()

! ==============================================================
! =======   Map the Extended mesh to the Regular grid   ========
! ==============================================================
   do i = 1,3
      dist = sqrt ( elvec(i,1)**2 + elvec(i,2)**2 + elvec(i,3)**2 )   ! get the minimal distance
      ipiv(i) = Rc_max / dist + 1  ! calc the trashold
      np(i) = ngrid(i) + 2*ipiv(i)    ! total points in extended mesh
   enddo ! do i
   emx1 = ipiv(1)   ! save offset of the extended mesh
   emx2 = ipiv(2)
   emx3 = ipiv(3)
   em1 = np(1)      ! save division of the extended mesh
   em2 = np(2)
   em3 = np(3)
   nem = em1*em2*em3
   write (*,*) '  ---------  Begin mesh Info --------- '
   write (*,*) ''
   write (*,'(a,f8.6)') ' Rc_max = ',Rc_max
   write (*,300) Ecut
   write (*,380) (g0(i),i=1,3)
   write (*,'(a,4i9)') ' Regular mesh: ',rm1, rm2, rm3, nrm
   write (*,'(a,f16.8)') 'dVol :',dvol
   write (*,'(a,4i9)') ' Extended mesh: ',em1, em2, em3, nem
   write (*,*) ' Elementary grid lattice vector :'
   write (*,2000)  (elvec(1,i),i=1,3)
   write (*,2000)  (elvec(2,i),i=1,3)
   write (*,2000)  (elvec(3,i),i=1,3)
   write (*,*) '  ---------   End mesh Info --------- '

! ==============================================================
! =======   Map the Extended mesh to the Regular grid   ========
! ==============================================================
   allocate (e2r(nem))
   do k = 0, em3-1    ! Loops over each axis
      do j = 0, em2-1
         do i = 0, em1-1
            ii = i - emx1    ! rest the offset
            jj = j - emx2
            kk = k - emx3
            ii = mod( ii + 100*rm1, rm1)            ! find the point in the regular mesh
            jj = mod( jj + 100*rm2, rm2)
            kk = mod( kk + 100*rm3, rm3)
            index = 1 + i + em1*j + em1*em2*k       ! calc index of the extended mesh point
            index0 = 1 + ii + rm1*jj + rm1*rm2*kk   ! calc index of the regular mesh point
            e2r(index) = index0                     ! save the value of the point in the regular mesh context
         enddo ! do i
      enddo ! do j
   enddo ! do k

! ==============================================================
! =======             Setup Atomic Mesh                 ========
! ==============================================================
! generate atomic mesh (sphere) within Rc_max
   write (*,*) 'atomic mesh'
   nam = 0
   do k = -emx3,emx3
      do j = -emx2,emx2
         do i = -emx1,emx1
            cvec(:) = elvec(1,:)*i + elvec(2,:)*j + elvec(3,:)*k
            dist = sqrt (cvec(1)**2 + cvec(2)**2 + cvec(3)**2)
            if ((Rc_max + drmax) .gt. dist) then   ! distance from point to mesh cell
               nam = nam + 1
            endif
         enddo ! do i
      enddo ! do j
   enddo ! do k

   allocate (am2rc(nam))       ! allocate arrays related to atomic mesh
   allocate (ram2rc(3,nam))
   allocate (ratom2g(3,natoms))
   index = 0
   do k = -emx3,emx3
      do j = -emx2,emx2
         do i = -emx1,emx1
            cvec(:) = elvec(1,:)*i + elvec(2,:)*j + elvec(3,:)*k
            dist = sqrt (cvec(1)**2 + cvec(2)**2 + cvec(3)**2)   ! distance to the mesh cell
            if ((Rc_max + drmax) .gt. dist) then                 ! ?? is the point within the Rc_max ??
               index = index + 1
               am2rc(index) = i + em1*j + em1*em2*k
               ram2rc(:,index) = cvec(:)
            endif
         enddo ! do i
      enddo ! do j
   enddo ! do k

   write (*,'(a,4i9)') ' Atomic mesh: ',2*emx1 + 1,2*emx2 + 1, 2*emx3 + 1, nam

! ==============================================================
! =======              SET UP FDM GRID                  ========
! ==============================================================
! set extended mesh (used for finite difference method)
   noff = 1
   mfd1 = rm1 + 2*noff
   mfd2 = rm2 + 2*noff
   mfd3 = rm3 + 2*noff
   nmfd = mfd1 * mfd2 * mfd3

! write out information about grids
   write (*,*) '  ----    the FDM grid informations   ----- '
   write (*,300) Ecut
   write (*,*) '  Normal FDM mesh: ',nrm
   write (*,400) rm1,rm2,rm3
   write (*,*) '  Extended FDM mesh: ', nmfd
   write (*,400) mfd1,mfd2,mfd3

! modified by honza
! set up auxilliary array keeping track of neighbors
! list of neighbors - 7,8 are currently redundant
   nneighij = 8
! relative coordinates
   neighij(1,1) = 1
   neighij(1,2) = -1
   neighij(2,1) = mfd1
   neighij(2,2) = -1*mfd1
   neighij(3,1) = mfd1*mfd2
   neighij(3,2) = -1*mfd1*mfd2
   neighij(4,1) = 0
   neighij(4,2) = 0

! recalculate the real space step of the mesh & convert it into a.u.
! remeber: only rectangular mesh can be used !!
! experimenting with (possibly) non-rectangular grid

! Vectors elvec(i,:) form the covariant base of the geometry.
! Vectors ervec(i,:) form the contravariant base of the geometry.
!
! dvol is actually kind of volume form.
! Explicitly: dvol = sqrt(det(metric)) (where metric would be the gramm-matrix (<elvec(i,:),elvec(j,:)>;i,j=1,3)
!

! 1st derivative coefficients
! No need to introduce any normalization factor - it's already present in the base vectors.
   d1f(1,1) = 1.0d0 / 2.0d0
   d1f(1,2) = -1.0d0 / 2.0d0
   d1f(2,1) = 1.0d0 / 2.0d0
   d1f(2,2) = -1.0d0 / 2.0d0
   d1f(3,1) = 1.0d0 / 2.0d0
   d1f(3,2) = -1.0d0 / 2.0d0
   d1f(4,1) = 0.0d0
   d1f(4,2) = 0.0d0

! end-modified by honza

! allocate arrays
   allocate (e2n (0:nmfd-1))
   allocate (n2e (0:nrm-1))

! allocate arrays
   allocate (vnaG (0:nrm-1))
   allocate (drhoG (0:nrm-1))
   allocate (rhoG0 (0:nrm-1))
   allocate (vcaG (0:nrm-1))
   allocate (vxcG (0:nrm-1))

   vnaG = 0.0d0
   drhoG = 0.0d0
   rhoG0 = 0.0d0
   vcaG = 0.0d0
   vxcG = 0.0d0

! --------------------------------------------------------------
!         set up mapping NORMAL(REGULAR) mesh INTO EXTENDED-FD mesh
! --------------------------------------------------------------
   index = 0
   index0 = mfd1*(mfd2+1)
   do k = 0, rm3-1
      do j = 0, rm2-1
         do i = 0, rm1-1
            index0 = index0 + 1
            n2e(index) = index0
!            write (*,*) ' norm : ',index,'ext : ',n2e(index)
            index = index + 1
         enddo ! do i
         index0 = index0 + 2*noff
      enddo ! do j
      index0 = index0 + 2*mfd1
   enddo ! do k

! --------------------------------------------------------------
!         Mapping EXTENDED-FD mesh TO NORMAL(REGULAR) mesh
! --------------------------------------------------------------
   do k = 0, mfd3-1
      do j = 0, mfd2-1
         do i = 0, mfd1-1
            ii = i - noff   ! rest the offset
            jj = j - noff
            kk = k - noff
            ii = mod( ii + 1000*rm1, rm1)   ! find the point in the basic mesh
            jj = mod( jj + 1000*rm2, rm2)
            kk = mod( kk + 1000*rm3, rm3)
            index =  i + mfd1*j + mfd1*mfd2*k    ! calc index of the fdm mesh point
            index0 = ii + rm1*jj + rm1*rm2*kk    ! calc index of the regular mesh point
            e2n(index) = index0                   ! save the value of the point in the basic mesh context
         enddo ! do i
      enddo ! do j
   enddo ! do k

! Format Statements
! ===========================================================================
100 format (2x, i5, 2x, a40, 2x, i2)
200 format (2x, 70('='))
225 format (2x,  3f15.7)
250 format (2x, ' Vol  = ', f16.6, ' [Ang^3] ')
300 format (2x, ' Ecut = ', f16.6, ' [Ry] ')
350 format (2x, ' dr = ', 3f16.8, ' [Ang]')
370 format (2x, ' dV = ', f16.8, ' [Ang^3]')
380 format (2x, ' Virtual initial grid point :', 3f16.8)
400 format (2x, ' Number of points = '3i9)
500 format (2x, ' Virtual center of mass : '3f16.8)
2000 format (2x, '  '3f16.8)

   return
end subroutine
