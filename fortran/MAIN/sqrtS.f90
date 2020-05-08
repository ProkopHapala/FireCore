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
subroutine sqrtS( Smat, norbitals, divide )
    use  workmat
    implicit none

! Arguments
    complex, dimension (:,:), pointer, intent(inout) :: Smat
    integer, intent (in) :: norbitals
    logical, intent(inout) :: divide  ! do we divide and conquer?

! globals
    real sqlami
    integer mineig
    integer norbitals_new
    integer info
    integer imu,jmu

    real*8, parameter :: overtol = 1.0d-4

! Procedure
! ===========================================================================
! Initialize some things

    if (divide) then
        call cheevd('V', 'U', norbitals, Smat, norbitals, slam, work, lwork, rwork , lrwork, iwork, liwork, info ) 
    else
        call cheev ('V', 'U', norbitals, Smat, norbitals, slam, xxxx, lwork, rwork , info)
    end if
    if (info .ne. 0) call diag_error (info, 0)

    ! xxxx = unused
    ! zzzz = Overlap eigenvectors
    ! yyyy = Hamiltonian

    ! CHECK THE LINEAR DEPENDENCE
    ! ****************************************************************************
    ! Fix the linear dependence problem. References: Szabo and Ostlund, Modern
    ! Quantum Chem. McGraw Hill 1989 p. 142; Szabo and Ostlund, Modern Quantum
    ! Chem. Dover 1996 p. 145. A tolerance for a small overlap eigenvalue is
    ! set by overtol.
    ! Determine the smallest active eigenvector
    mineig = 0
    do imu = 1, norbitals
        if (slam(imu) .lt. overtol) mineig = imu
    end do

    ! You can specify a specific number of orbitals to drop with this next line, by uncommenting it.
    ! mineig = 0  {Don't drop any}
    mineig = mineig + 1
    norbitals_new = norbitals + 1 - mineig
    if (norbitals_new .ne. norbitals) then
        write (*,*) '  '
        write (*,*) ' WARNING. ############################ '
        write (*,*) ' Linear dependence encountered in basis set. '
        write (*,*) ' An overlap eigenvalue is very small. '
        write (*,*) norbitals - norbitals_new, ' vectors removed. '
        write (*,*) ' Spurious orbital energies near zero will '
        write (*,*) ' appear as a result of dropping these orbitals'
        write (*,*) ' You can change this by adjusting overtol in '
        write (*,*) ' kspace.f '
        write (*,*) '  '
        ! if(ishort .eq. 1) then    ! Don't print out again if done above
        !     write (*,*) '            The overlap eigenvalues: '
        !     write (*,*) ' ********************************************** '
        !     write (*,200) (slam(imu), imu = 1, norbitals)
        !     else                      ! They asked for extra printout
        !     write(*,*) ' '
        !     write(*,*) ' Eigenvectors that correspond to eigenvalues'
        !     write(*,*) ' that were eliminated.  These might provide'
        !     write(*,*) ' insight into what atomic orbitals are causing'
        !     write(*,*) ' the problem.'
        !     write(*,*) ' '
        !     do imu = 1, mineig - 1
        !         write(*,*) ' eigenvector',imu
        !         do jmu = 1, norbitals
        !             write(*,*) jmu,' ',zzzz(jmu,imu)
        !         end do
        !     end do
        ! end if
        write (*,*) ' '
        do imu = mineig, norbitals
            jmu = imu - mineig + 1
            Smat(:,jmu) = Smat(:,imu)
            slam(jmu) = slam(imu)
        end do
    end if ! (norbitals_new .ne. norbitals)

    ! CALCULATE (S^-1/2) --> sm1
    ! ****************************************************************************
    ! In a diagonal reperesentation (Udagger*S*U = s, s is a diagonal matrix)
    ! We just take the inverse of the square roots of the eigenvalues to get
    ! s^-1/2. Then we 'undiagonalize' the s^-1/2 matrix back to get
    ! S^-1/2 = U*s^-1/2*Udagger.
    ! Note: We do S^-1/4 here, because the sqlami contribution get squared
    ! after it is combined with overlap.
    do imu = 1, norbitals_new
        sqlami = slam(imu)**(-0.25d0)
        Smat(:,imu) = Smat(:,imu)*sqlami
    end do

    call cgemm ('N', 'C', norbitals, norbitals, norbitals_new, cmplx(1.0d0,0.0d0), Smat, norbitals, Smat, norbitals, cmplx(0.0d0,0.0d0), xxxx, norbitals)

end subroutine sqrtS

