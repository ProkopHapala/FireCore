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
    use debug
    use  workmat
    use interactions, only: norbitals_new
    use options, only: idebugWrite
    implicit none

! Arguments
    complex, dimension (norbitals,norbitals), intent(inout) :: Smat
    integer, intent(in) :: norbitals
    logical, intent(in) :: divide  ! do we divide and conquer?

! globals
    real sqlami
    integer mineig
    !integer norbitals_new
    integer info
    integer imu,jmu

    integer ifile

    real*8, parameter :: overtol = 1.0d-4

    complex, dimension (1) :: work_query

! Procedure
! ===========================================================================
! DIAGONALIZE THE OVERLAP MATRIX
! ****************************************************************************
! If you are within the scf loop, you do not have to recalculate the overlap.
! NPA

    ! Call the diagonalizer
    !         if (wrtout) then
    !           write (*,*) ' Call diagonalizer for overlap. '
    !           write (*,*) '                  The overlap eigenvalues: '
    !           write (*,*) ' ******************************************************* '
    !!         end if


    zzzz = Smat

    !write (*,*) "!!!! DEBUG sqrtS debug_writeMatFile(Sk_sqrtS.log) norbitals=",norbitals, ' lwork = ',lwork, ' lrwork = ',lrwork, ' liwork = ',liwork, ' divide = ',divide
    !call debug_writeMatFile_cmp( "Sk_sqrtS_", zzzz, norbitals, norbitals, 0 )

    if (divide) then
        call zheevd('V', 'U', norbitals, zzzz, norbitals, slam, work,  lwork, rwork , lrwork, iwork, liwork, info )
    else
! first find optimal working space
        !  zheev (JOBZ,UPLO,N,A[lda,*],lda,W_evals,WORK,LWORK,RWORK,integer 	INFO )		
        call zheev ('V', 'U', norbitals, zzzz, norbitals, slam, work_query,   -1, rwork, info)
        if (info .ne. 0) call diag_error (info, 0)
        !call debug_writeMatFile_cmp( "S_query_", zzzz, norbitals, norbitals, 0 )

! resize working space
        lwork = work_query(1)
        deallocate (work)
        allocate   (work(lwork))
! diagonalize the overlap matrix with the new working space
        call zheev ('V', 'U', norbitals, zzzz, norbitals, slam, work,  lwork, rwork , info)
    end if

    if (info .ne. 0) call diag_error (info, 0)
!    if (ishort .eq. 1 .and. wrtout) then
!          write (*,100) slam(1), slam(norbitals)
!    else if (wrtout) then
!          write (*,200) (slam(imu), imu = 1, norbitals)
!    end if

! xxxx = unused
! zzzz = Overlap eigenvectors
! yyyy = Hamiltonian

    !write(*,*) "!!!! DEBUG S-eigenvalues: ", slam(:)

    !write (*,*) "!!!! DEBUG sqrtS() debug_writeMatFile(zzzz_pre.log) norbitals=",norbitals, " norbitals_new= ", norbitals_new, "info ", info
    !call debug_writeMatFile_cmp( "zzzz_pre1_", zzzz, norbitals, norbitals, 0 )

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

! You can specify a specific number of orbitals to drop with this
! next line, by uncommenting it.
! mineig = 0  {Don't drop any}

    mineig = mineig + 1
    norbitals_new = norbitals + 1 - mineig

    if ( (norbitals_new .ne. norbitals) .or.  (idebugWrite .gt. 0)  ) then
        write (*,*) '  '
        write (*,*) ' !!!!!!!! WARNING. !!!!!!!!! '
        write (*,*) ' norbitals_new(',norbitals_new,') .ne. norbitals(',norbitals,')'
        write (*,*) ' => Linear dependence encountered in basis set. An overlap eigenvalue is very small.', norbitals - norbitals_new, ' vectors removed. '
        write (*,*) ' Spurious orbital energies near zero will appear as a result of dropping these orbitals'
        write (*,*) ' You can change this by adjusting overtol(',overtol,') in kspace.f '
        do imu = 1, norbitals
            write (*,*) 'slam[',imu,']=',slam(imu),' vs overtol(',overtol,')'
        end do

        !if(ishort .eq. 1) then    ! Don't print out again if done above
        !    write (*,*) '            The overlap eigenvalues: '
        !    write (*,*) ' ********************************************** '
        !    write (*,200) (slam(imu), imu = 1, norbitals)
        !else                      ! They asked for extra printout
        !    write(*,*) ' '
        !    write(*,*) ' Eigenvectors that correspond to eigenvalues'
        !    write(*,*) ' that were eliminated.  These might provide'
        !    write(*,*) ' insight into what atomic orbitals are causing'
        !    write(*,*) ' the problem.'
        !    write(*,*) ' '
        !    do imu = 1, mineig - 1
        !    write(*,*) ' eigenvector',imu
        !    do jmu = 1, norbitals
        !        write(*,*) jmu,' ',zzzz(jmu,imu)
        !    end do
        !    end do
        !end if ! ishort
        write (*,*) ' '

        do imu = mineig, norbitals
            jmu = imu - mineig + 1
            zzzz(:,jmu) = zzzz(:,imu)
            slam(jmu) = slam(imu)
        end do
    end if ! norbitals_new .ne. norbitals

    !write (*,*) "!!!! DEBUG sqrtS() debug_writeMatFile(zzzz_pre.log) norbitals=",norbitals, " norbitals_new= ", norbitals_new
    !call debug_writeMatFile_cmp( "zzzz_pre2_", zzzz, norbitals, norbitals, 0 )

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
        zzzz(:,imu) = zzzz(:,imu)*sqlami
    end do

    !write (*,*) "!!!! DEBUG sqrtS() debug_writeMatFile(zzzz_pre.log) norbitals=",norbitals, " norbitals_new= ", norbitals_new
    !call debug_writeMatFile_cmp( "zzzz_pre3_", zzzz, norbitals, norbitals, 0 )

    call zgemm ('N', 'C', norbitals, norbitals, norbitals_new, cmplx(1.0d0,0.0d0), zzzz, norbitals, zzzz, norbitals, cmplx(0.0d0,0.0d0), xxxx, norbitals)

    !write (*,*) "!!!! DEBUG sqrtS() debug_writeMatFile(zzzz.log) norbitals=",norbitals
    !call debug_writeMatFile_cmp( "zzzz_", zzzz, norbitals, norbitals, 0 )
    !write (*,*) "!!!! DEBUG sqrtS() debug_writeMatFile(Sm12.log) norbitals=",norbitals
    !call debug_writeMatFile_cmp( "xxxx_", xxxx, norbitals, norbitals, 0 )

end subroutine sqrtS