! copyright info:
!
!                             @Copyright 2001
!                           Fireball Committee
! University of Utah - James P. Lewis, Chair
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
! Ohio State University - Dave Drabold
 
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

 
! allocate_rho.f90
! Program Description
! ===========================================================================
!       This routine allocates the arrays which store information about 
! wavefunctions, neutral atomic potentials (vna) of each specie and finally
! the arrays of numerical grid 
!
! ===========================================================================
! Code written by:

! ===========================================================================
!
! Program Declaration
! ===========================================================================
  subroutine allocate_grid !(natoms, nspecies)
    use configuration, only: natoms, nspecies
    use wavefunction
    use vnneutral
    use interactions
    use grid
    use options, only : verbosity
    implicit none
    logical :: need_realloc

! Argument Declaration and Description
! ===========================================================================
! Inputinitcharges_ks_
!    integer, intent (in) :: natoms
!    integer, intent (in) :: nspecies
 
! Local Parameters and Data Declaration
! ===========================================================================
 
! Local Variable Declaration and Description
! ===========================================================================
 
! Allocate Arrays
! ===========================================================================
 
! Procedure
! ===========================================================================
    !write(*,*) "subroutine reallocate_grid() verbosity",verbosity
    if(verbosity.gt.0) write(*,*) "subroutine reallocate_grid() "
    need_realloc = allocated(mesh_wf)
    if (need_realloc) then
        if (verbosity.gt.0) write(*,*) "WARNING allocate_grid(): freeing existing grid buffers"
        call deallocate_grid()
    end if

! wave function part
    allocate (mesh_wf (nsh_max, nspecies))
    allocate (drr_wf (nsh_max, nspecies))
    allocate (rmax_wf (nsh_max, nspecies))
    allocate (rr_wf (wfmax_points, nsh_max, nspecies))
    allocate (wf_spline (wfmax_points, nsh_max, nspecies))
    allocate (wf (wfmax_points, nsh_max, nspecies))

! vneutral part
    allocate (mesh_na (nspecies))
    allocate (rmax_na (nspecies))
    allocate (drr_na (nspecies))
    allocate (rr_na (max_vna_points, nspecies))
    allocate (vnna_spline (max_vna_points, nspecies))
    allocate (vnna (max_vna_points, nspecies))

! Deallocate Arrays
! ===========================================================================
 
! Format Statements
! ===========================================================================
 
    return
  end subroutine allocate_grid

  subroutine deallocate_grid()
    use configuration, only: natoms, nspecies
    use wavefunction
    use vnneutral
    use interactions
    use grid
    use options, only : verbosity
    implicit none

    if(verbosity.gt.0) write(*,*) "deallocate_grid()"

! wave function part
    if (allocated(mesh_wf))       deallocate(mesh_wf)
    if (allocated(drr_wf))        deallocate(drr_wf)
    if (allocated(rmax_wf))       deallocate(rmax_wf)
    if (allocated(rr_wf))         deallocate(rr_wf)
    if (allocated(wf_spline))     deallocate(wf_spline)
    if (allocated(wf))            deallocate(wf)

! vneutral part
    if (allocated(mesh_na))       deallocate(mesh_na)
    if (allocated(rmax_na))       deallocate(rmax_na)
    if (allocated(drr_na))        deallocate(drr_na)
    if (allocated(rr_na))         deallocate(rr_na)
    if (allocated(vnna_spline))   deallocate(vnna_spline)
    if (allocated(vnna))          deallocate(vnna)

  end subroutine deallocate_grid
