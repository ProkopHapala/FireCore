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
subroutine ktransform( kpoint, n, Sk, Hk )
        use configuration  
        use density
        use dimensions
        use interactions
        use neighbor_map
        implicit none

! Argument Declaration and Description
! ===========================================================================
! Input
!        integer, intent (in) :: icluster
        !integer, intent (in) :: ikpoint
         real,    dimension (3), intent (in) :: kpoint
         integer,                intent (in) :: n
         complex, dimension (n,n), intent (out) :: Sk
         complex, dimension (n,n), intent (out) :: Hk
!        integer, intent (in) :: iqout
!        integer, intent (in) :: iwrteigen
!        integer, intent (in) :: Kscf
!        integer, intent (in) :: nprocs
!        integer, intent (in) :: my_proc
!        integer, intent (in) :: nkpoints
!        integer, intent (in) :: iwrtdos
!        integer, intent (in) :: iwrthop
!        integer, intent (in) :: iwrtatom
!        integer, intent (in) :: itrans

! Local Variable Declaration and Description
! ===========================================================================
        integer iatom
        integer imu
        integer inu
        integer in1, in2
        integer ineigh
        integer jatom
        integer jmu
        integer jnu
        integer mbeta

        real dot

        real, dimension (3) :: vec

        complex phase

! ====== Body


! We built H/S(K) going through the list of atoms and their neighbors.
! How do we know where the particular (i,j) spot belongs in the grand
! matrix?  This should help you understand the degelec shifting business:
!
!                  atom 1         atom2            atom  3
!
!              nu1 nu2 nu3   nu1 nu2 nu3 nu4     nu1 nu2 nu3
!
!                           _________________
!         mu1               |               |
!  atom 1 mu2               |    H(1,2)     |
!         mu3               |               |
!                           -----------------
!         mu1
!  atom 2 mu2
!         mu3
!         mu4
!
!         mu1
!  atom3  mu2
!         mu3
!
! to the (1,2) portion at the right place we use degelec(iatom), which is
! passed, it remembers how many orbitals one must skip to get to the spot
! reserved for iatom, e.g. in this example degelec(1)=0, degelc(2)=3.

!
! COMPUTE S(K) AND H(K)
! ****************************************************************************
! Find the overlap and Hamiltonian matrices s(mu,nu,i,j) and h(mu,nu,i,j)
! in k space. Here iatom is an atom in the central cell and jatom a
! neighboring atom. The k-space matrix is found from the real space matrix by:
! s(k) = sum(l) exp(iatom*k - dot - (l + bj - bi)) s((0,iatom), (l,jatom))
! h(k) = sum(l) exp(iatom*k - dot - (l + bj - bi)) h((0,iatom), (l,jatom))
! Initialize to zero first

        Sk = complex(0.0,0.0)
        Hk = complex(0.0,0.0)
        ! We now loop over all neighbors jatom of iatom.
        do iatom = 1, natoms
         in1 = imass(iatom)
         
         ! --- H_mat and S_mat pair interactions
         do ineigh = 1, neighn(iatom)
          mbeta = neigh_b(ineigh,iatom)
          jatom = neigh_j(ineigh,iatom)
          in2 = imass(jatom)
          ! So this matrix element goes in the i,j slot.
          vec(:) = xl(:,mbeta) + ratom(:,jatom) - ratom(:,iatom)
          dot = kpoint(1)*vec(1) + kpoint(2)*vec(2) + kpoint(3)*vec(3)
          phase = cmplx(cos(dot),sin(dot))
          do inu = 1, num_orb(in2)
           jnu = inu + degelec(jatom)
           do imu = 1, num_orb(in1)
            jmu = imu + degelec(iatom)
            Sk(jmu,jnu) = Sk(jmu,jnu) + phase*s_mat(imu,inu,ineigh,iatom)   ! Sk === zzzz
            Hk(jmu,jnu) = Hk(jmu,jnu) + phase*h_mat(imu,inu,ineigh,iatom)   ! Hk === yyyy
           end do ! do inu
          end do ! do imu
         end do ! do ineigh

         ! --- VNL pair interactions
         do ineigh = 1, neighPPn(iatom)
          mbeta = neighPP_b(ineigh,iatom)
          jatom = neighPP_j(ineigh,iatom)
          in2 = imass(jatom)
          ! So this matrix element goes in the i,j slot.
          vec(:) = xl(:,mbeta) + ratom(:,jatom) - ratom(:,iatom)
          dot = kpoint(1)*vec(1) + kpoint(2)*vec(2) + kpoint(3)*vec(3)
          phase = cmplx(cos(dot),sin(dot))
          do inu = 1, num_orb(in2)
           jnu = inu + degelec(jatom)
           do imu = 1, num_orb(in1)
            jmu = imu + degelec(iatom)
            Hk(jmu,jnu) = Hk(jmu,jnu) + phase*vnl(imu,inu,ineigh,iatom)
           end do ! do imu
          end do ! do inu

         end do ! do inegh
        end do ! do iatom

        return
      end subroutine ktransform

