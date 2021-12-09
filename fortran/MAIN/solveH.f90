! copyright info:
!
!                             @Copyright 2001
!                           Fireball Committee
! Brigham Young University - James P. Lewis, Chair
! Arizona State University - Otto F. Sankey
! University of Regensburg - Juergen Fritsch
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


! kspace_blas.f90
! Program Description
! ===========================================================================
!       This is a version of kspace.f that uses the blas library
! ===========================================================================
! Original code written by Otto F. Sankey with modification by Alex A. Demkov
! and Jose Ortega

! Code rewritten by:
! James P. Lewis
! Kurt R. Glaesemann
! Henry Eyring Center for Theoretical Chemistry
! Department of Chemistry
! University of Utah
! 315 S. 1400 E.
! Salt Lake City, UT 84112-0850
! FAX 801-581-4353
! Office telephone 801-585-1078
! ===========================================================================
!
! Program Declaration
! ===========================================================================
!subroutine kspace (nprocs, my_proc, Kscf, iqout, icluster, iwrteigen, ikpoint, sks, nkpoints, iwrtdos, iwrthop, iwrtatom, itrans)

subroutine solveH( ikpoint, kpoint )
        use loops
        use kpoints
        use options
        use configuration  
        use density
        use dimensions
        use interactions
        use neighbor_map
        use workmat
        use debug
        implicit none

! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: ikpoint
        real,  dimension (3), intent (in) :: kpoint

! Local Parameters and Data Declaration
! ===========================================================================
        real, parameter :: overtol = 1.0d-4

! Local Variable Declaration and Description
! ===========================================================================
        !integer iatom
        integer imu, inu, jmu
        integer info
        !integer inu
        !integer in1, in2
        !integer ineigh
        !integer ishort
        !integer jatom
        !integer jmu
        !integer jnu
        !integer mbeta

        real dot

        real, dimension (norbitals) :: eigen
        !real, dimension (norbitals) :: slam
        !real, dimension (3) :: vec

        complex a0
        complex a1
        complex phase

! A bunch of memory to be used in many ways
!        complex, dimension (:, :), allocatable :: xxxx
        complex, dimension (:, :), allocatable :: Hk
        complex, dimension (:, :), allocatable :: Sk
!        complex, dimension (:, :, :), allocatable, save :: sm12_save

! work vector for cheev/cheevd
!        complex, allocatable, dimension (:) :: work
!        real,    allocatable, dimension (:) :: rwork
!        integer, allocatable, dimension (:) :: iwork
!        integer lwork, lrwork, liwork
        logical, parameter :: divide = .false.  ! do we divide and conquer?

! Procedure
! ===========================================================================
! Initialize some things

        if(idebugWrite .gt. 0) write(*,*) "DEBUG solveH()"

        !ishort = 1
        !if (iwrteigen .eq. 1) ishort = 0

        allocate ( xxxx(norbitals,norbitals) )
        allocate ( Hk(norbitals,norbitals) )
        allocate ( Sk(norbitals,norbitals) )
        allocate ( slam(norbitals) )

        call debug_writeBlockedMat( "Sk.log", Sk )
        call debug_writeBlockedMat( "Hk.log", Hk )
        stop  ! DEBUG

        call ktransform( kpoint, norbitals, Sk, Hk )
        if (divide) then
            lwork  = 100*norbitals + norbitals*norbitals
            lrwork = 100*norbitals + 3*norbitals*norbitals ! for old versions of cheevd
            liwork = 10*norbitals
            allocate (work(lwork))
            allocate (iwork(liwork))
            allocate (rwork(lrwork))
        else
            lwork = norbitals*norbitals ! Use xxxx, Hk and Sk for work area
            lrwork = 3*norbitals
            allocate (rwork(lrwork))
        end if
        if (.not. allocated(sm12_save)) then
            allocate (sm12_save(norbitals,norbitals,nkpoints))
        end if
        if (Kscf .eq. 1) then
            call sqrtS( Sk, norbitals, divide )
            do inu = 1, norbitals
                do imu = 1, norbitals
                    sm12_save(imu,inu,ikpoint) = xxxx(imu,inu)
                end do
            end do
        else    ! Kscf .eq. 1 
            xxxx(:,:) = sm12_save(:,:,ikpoint)  ! Restore S^-1/2 from s(k)^-1/2,
        end if  ! Kscf .eq. 1 

! xxxx = S^-1/2 in AO basis
! Sk = Unused (used as temporary work area below)
! Hk = Hamiltonian in AO basis

! CALCULATE (S^-1/2)*H*(S^-1/2)
! ****************************************************************************
        call chemm ( 'R', 'U', norbitals, norbitals, a1, xxxx, norbitals, Hk, norbitals, a0, Sk, norbitals ) ! Set M=H*(S^-.5)
        call chemm ( 'L', 'U', norbitals, norbitals, a1, xxxx, norbitals, Sk, norbitals, a0, Hk, norbitals ) ! Set Z=(S^-.5)*M
! xxxx = S^-1/2 in AO basis
! Sk = Unused (used as complex workspace in cheev call below)
! Hk = Hamiltonian in the MO basis set


! DIAGONALIZE THE HAMILTONIAN.
! ****************************************************************************

! Eigenvectors are needed to calculate the charges and for forces!
        if (divide) then
          call cheevd('V', 'U', norbitals, Hk, norbitals, eigen, work, lwork, rwork , lrwork, iwork, liwork, info ) 
        else
          call cheev ('V', 'U', norbitals, Hk, norbitals, eigen, Sk, lwork, rwork, info)
        end if
        if (info .ne. 0) call diag_error (info, 0)

! FIXME - Should only go up to norbitals_new, but we do not know what
! eigenvalues and eigenvectors correspond to the removed MO's.  Their
! eigenvalues will be very close to zero, but not exactly.  Also, we do not
! know if a real eigen value is near zero.

! INFORMATION FOR THE LOWDIN CHARGES
! ****************************************************************************
! xxxx = S^-1/2 in AO basis
! Sk = Unused
! Hk = Hamiltonian eigenvectors in the MO basis
        eigen_k(1:norbitals,ikpoint) = eigen(:)
        if (iqout .ne. 2) blowre(:,:,ikpoint) = real(Hk(:,:))
        if (iqout .ne. 2 .and. icluster .ne. 1) blowim(:,:,ikpoint) = aimag(Hk(:,:))
        call chemm ( 'L', 'U', norbitals, norbitals, a1, xxxx, norbitals, Hk, norbitals, a0, Sk, norbitals )
        bbnkre(:,:,ikpoint) = real(Sk(:,:))
        if (icluster .ne. 1) bbnkim(:,:,ikpoint) = aimag(Sk(:,:))

! We did a symmetric orthogonalization followed by a diagnalization
! of the Hamiltonian in this "MO" basis set. This yields a net
! canonical diagnolization with matrix bbnk.

! xxxx = S^-1/2 in AO basis
! Sk = S^-1/2 * Hk
! Hk = Hamiltonian eigenvectors in the MO basis

! Deallocate Arrays
! ===========================================================================
        deallocate (xxxx)
        deallocate (Hk)
        deallocate (Sk)
        deallocate (slam)
        deallocate (rwork)
        if (divide) then
          deallocate (work)
          deallocate (iwork)
        end if

! Format Statements
! ===========================================================================
100     format (2x, ' eigenvalue(1) = ', f10.6, ' eigenvalue(norbitals) = ', f10.6)
200     format (4(2x, f12.4))
        return
      end subroutine solveH

! ============ diag_error 
subroutine diag_error (info, istyle)
        implicit none
    ! Arguments
        integer, intent (in) :: info
        integer, intent (in) :: istyle
    ! Body
        write (*,*) '  '
        write (*,*) ' Diagonalization not successful, info = ', info
        if (info .lt. 0) then
            write (*,*) ' The ', info, '-th argument had an illegal '
            write (*,*) ' value. '
        else if(istyle .eq. 0)then
    ! LAPACK style errors
            write (*,*) ' It failed to converge.'
            write (*,*) info, ' off-diagonal elements of an intermediate'
            write (*,*) ' tridiagonal form did not converge to zero. '
    ! SCALAPACK style errors
        else if (mod(info,2) .ne. 0) then
            write (*,*) ' one or more eigenvectors failed to converge.'
            write (*,*) ' This should not have occured.  Send e-mail to'
            write (*,*) ' scalapack@cs.utk.edu if you feel like it'
        else if (mod(info/2,2) .ne. 0) then
            write (*,*) ' DARNGER Will Robinson.  Mr. Smith is in the house. '
            write (*,*) ' eigenvectors corresponding to one or more clusters '
            write (*,*) ' of eigenvalues could not be reorthogonalized'
            write (*,*) ' because of insufficient workspace.'
            write (*,*) ' We will blindly go on and hope for the best. '
            return
        else if (mod(info/4,2) .ne. 0) then
            write (*,*) ' space limit prevented computing all of the eigenvectors '
        else if (mod(info/8,2) .ne. 0) then
            write (*,*) ' PCSTEBZ failed to compute eigenvalues '
            write (*,*) ' This is very strange indeed.  It should not happen'
            write (*,*) ' Send e-mail to scalapack@cs.utk.edu if you feel like it'
        end if
        stop
    end subroutine diag_error
