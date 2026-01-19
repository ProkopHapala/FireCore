







module debug
    implicit none


! --------------------------
! DEBUG : TO EXPORT For checking /pyBall/FireballOCL/OCL_Hamiltonian.py
! AvgRho diagnostics (filled in average_rho.f90 when enabled)
         logical :: diag_avg_rho_enable = .false.
         integer :: diag_avg_rho_iatom = -1
         integer :: diag_avg_rho_jatom_target = -1
         integer :: diag_avg_rho_mbeta_target = -1
         integer :: diag_avg_rho_ineigh = -1
         integer :: diag_avg_rho_in1 = -1
         integer :: diag_avg_rho_in2 = -1
         integer :: diag_avg_rho_jatom = -1
         integer :: diag_avg_rho_mbeta = -1
         real, dimension(3,3) :: diag_avg_rho_eps2c = 0.0
         real, dimension(6,6) :: diag_avg_rho_sm = 0.0
         real, dimension(6,6) :: diag_avg_rho_rhom2c = 0.0
         real, dimension(6,6) :: diag_avg_rho_rhom3c = 0.0
         real, dimension(8,8) :: diag_avg_rho_rhooff_3c = 0.0
         real, dimension(8,8) :: diag_avg_rho_rhooff_final = 0.0

! VCA component debug exports (filled in assemble_ca_2c.f90 when idebugWrite>0)
! These arrays are NOT part of the production state (kept in MODULES/debug.f90 only).
! They exist only to provide Fortran reference component breakdown for PyOpenCL parity tests.
         real, dimension (:, :, :, :), allocatable :: vca_ontopl
         real, dimension (:, :, :, :), allocatable :: vca_ontopr
         real, dimension (:, :, :, :), allocatable :: vca_atom

! VCA trace selector (consumed by assemble_ca_2c.f90 when idebugWrite>0)
! Allows targeting a specific directed neighbor (iatom,jatom,mbeta) and shell (isorp).
! Use -1 as wildcard for any field.
         logical :: diag_vca_enable = .false.
         integer :: diag_vca_iatom  = -1
         integer :: diag_vca_jatom  = -1
         integer :: diag_vca_mbeta  = -1
         integer :: diag_vca_isorp  = -1

! EWALDSR component debug exports (filled in assemble_ca_2c.f90 / assemble_ca_3c.f90 when idebugWrite>0)
! These arrays are NOT part of the production state (kept in MODULES/debug.f90 only).
! They exist only to provide Fortran reference component breakdown for PyOpenCL parity tests.
         real, dimension (:, :, :, :), allocatable :: ewaldsr_2c_atom
         real, dimension (:, :, :, :), allocatable :: ewaldsr_2c_ontop
         real, dimension (:, :, :, :), allocatable :: ewaldsr_3c

! EWALDSR trace selector (consumed by assemblers when idebugWrite>0)
! Allows targeting a specific directed neighbor (iatom,jatom,mbeta) and third center (ialp).
! Use -1 as wildcard for any field.
         logical :: diag_ewald_enable = .false.
         integer :: diag_ewald_iatom  = -1
         integer :: diag_ewald_jatom  = -1
         integer :: diag_ewald_mbeta  = -1
         integer :: diag_ewald_ialp   = -1

! --------------------------






    contains


subroutine debug_vca_prepare()
    use configuration, only: natoms
    use neighbor_map,  only: neigh_max
    use interactions,  only: numorb_max
    use options,       only: idebugWrite, verbosity
    implicit none
    if (idebugWrite .le. 0) return
    if (.not. allocated(vca_ontopl)) allocate( vca_ontopl(numorb_max, numorb_max, neigh_max, natoms) )
    if (.not. allocated(vca_ontopr)) allocate( vca_ontopr(numorb_max, numorb_max, neigh_max, natoms) )
    if (.not. allocated(vca_atom  )) allocate( vca_atom  (numorb_max, numorb_max, neigh_max, natoms) )
    vca_ontopl = 0.0d0
    vca_ontopr = 0.0d0
    vca_atom   = 0.0d0
    if (verbosity .gt. 0) write(*,*) 'DEBUG VCA export enabled: allocated+zeroed vca_ontopl/vca_ontopr/vca_atom'
end subroutine debug_vca_prepare


subroutine debug_vca_zero()
    use options, only: idebugWrite
    implicit none
    if (idebugWrite .le. 0) return
    if (allocated(vca_ontopl)) vca_ontopl = 0.0d0
    if (allocated(vca_ontopr)) vca_ontopr = 0.0d0
    if (allocated(vca_atom  )) vca_atom   = 0.0d0
end subroutine debug_vca_zero


subroutine debug_ewald_prepare()
    use configuration, only: natoms
    use neighbor_map,  only: neigh_max
    use interactions,  only: numorb_max
    use options,       only: idebugWrite, verbosity
    implicit none
    if (idebugWrite .le. 0) return
    if (.not. allocated(ewaldsr_2c_atom )) allocate( ewaldsr_2c_atom (numorb_max, numorb_max, neigh_max, natoms) )
    if (.not. allocated(ewaldsr_2c_ontop)) allocate( ewaldsr_2c_ontop(numorb_max, numorb_max, neigh_max, natoms) )
    if (.not. allocated(ewaldsr_3c     )) allocate( ewaldsr_3c      (numorb_max, numorb_max, neigh_max, natoms) )
    if (verbosity .gt. 0) write(*,*) 'DEBUG EWALDSR export enabled: allocated+zeroed ewaldsr_2c_atom/ewaldsr_2c_ontop/ewaldsr_3c'
end subroutine debug_ewald_prepare


subroutine debug_ewald_zero()
    use options, only: idebugWrite
    implicit none
    if (idebugWrite .le. 0) return
    if (allocated(ewaldsr_2c_atom )) ewaldsr_2c_atom  = 0.0d0
    if (allocated(ewaldsr_2c_ontop)) ewaldsr_2c_ontop = 0.0d0
    if (allocated(ewaldsr_3c     )) ewaldsr_3c       = 0.0d0
end subroutine debug_ewald_zero






subroutine debug_writeArray_1i( name, arr, n )
    ! ===== Parameters
    character(len=*)      , intent(in)   :: name
    integer, dimension( : ), intent (in) :: arr
    integer,                 intent (in) :: n
    ! ===== Variables
    integer i
    ! ===== Body
    do i=1,n
        write(*,*) name,"[",i,"]: ", arr(i)
    end do
end subroutine debug_writeArray_1i

subroutine debug_writeAtomCartes( ifile, name, arr, n )
    ! ===== Parameters
    integer,                  intent (in) :: ifile
    character(len=*)        , intent(in)  :: name
    real, dimension( :, : ) , intent (in) :: arr
    integer,                  intent (in) :: n
    ! ===== Variables
    integer i
    ! ===== Body
    do i=1,n
        if(ifile .eq. 0) then
            write(*,*) name,"[",i,"]: ", arr(1,i),arr(2,i),arr(3,i)
        else
            write(ifile,*) name,"[",i,"]: ", arr(1,i),arr(2,i),arr(3,i)
        end if
    end do
end subroutine debug_writeAtomCartes


subroutine debug_writeMat( ifile, mat, n, m )
    ! ===== Parameters
    integer,                 intent (in) :: ifile
    real, dimension( :, : ), intent (in) :: mat
    integer,                 intent (in) :: n, m
    ! ===== Variables
    integer i, j
    ! ===== Body
    do i=1,n
        do j=1,m
            write(ifile,'(e20.10)',advance='no' ) mat(i,j) 
        end do
        write (ifile,*)
    end do
end subroutine debug_writeMat

subroutine debug_writeMatFile( fname, mat, n, m )
    ! ===== Parameters
    character(len=*), intent(in) :: fname
    real, dimension( :, : ), intent (in) :: mat
    integer,                 intent (in) :: n, m
    ! ===== Variables
    integer                 ifile
    ! ===== Body
    ifile = 11111
    open( ifile, file=fname, status='unknown' )
    call debug_writeMat( ifile, mat, n, m )
    close(ifile)
end subroutine debug_writeMatFile

subroutine debug_writeMat_cmp( ifile, mat, n, m )
    ! ===== Parameters
    integer,                 intent (in) :: ifile
    complex, dimension( :, : ), intent (in) :: mat
    integer,                 intent (in) :: n, m
    ! ===== Variables
    integer i, j
    ! ===== Body
    do i=1,n
        do j=1,m
            write(ifile,'(e20.10,e20.10)',advance='no' ) real(mat(i,j)),aimag(mat(i,j)) 
        end do
        write (ifile,*)
    end do
end subroutine debug_writeMat_cmp

subroutine debug_writeMatFile_cmp( fname, mat, n, m, iname )
    ! ===== Parameters
    character(len=*), intent(in) :: fname
    complex, dimension( :, : ), intent (in) :: mat
    integer,                 intent (in) :: n, m, iname
    ! ===== Variables
    integer                 ifile
    character(40)   :: fname_
    character(4)    :: name
    ! ===== Body
    ifile = 11111
    write (name,'(i4.4)') iname
    fname_ = trim(fname)//trim(name)//".log"
    open( ifile, file=fname_, status='unknown' )
    call debug_writeMat_cmp( ifile, mat, n, m )
    close(ifile)
end subroutine debug_writeMatFile_cmp

subroutine debug_writeBlockedMat( name, M, unit  )
    use configuration
    use neighbor_map
    use interactions
    use, intrinsic :: iso_fortran_env
    ! allocate (s_mat (numorb_max, numorb_max, neigh_max, natoms))
    ! ===== Parameters
    character(len=*)                                               , intent(in)  :: name
    real, dimension( numorb_max, numorb_max, neigh_max, natoms ), intent (in) :: M
    ! ===== Variables
    integer iatom, jatom, ineigh, inu, imu,  jnu, jmu, in1,in2
    integer ioff, joff, nnu, nmu, nng
    integer unit
    logical is_file_output
    ! ===== Body

    write(*,*) "Standard output unit:", OUTPUT_UNIT

    is_file_output = (unit /= OUTPUT_UNIT)
    if (is_file_output) then
        open(unit, file=name, status='unknown')
    else
        write(*,*) "  "
        write(*,*) "debug_writeBlockedMat(name= ", name, " )"
    end if
    do iatom = 1, natoms
        in1  = imass(iatom)
        ioff = degelec(iatom)
        nmu  = num_orb(in1)
        nng  = neighn(iatom)
        do ineigh = 1, nng
         !mbeta = neigh_b(ineigh,iatom)
         jatom = neigh_j(ineigh,iatom)
         in2   = imass(jatom)
         joff  = degelec(jatom)
         nnu   = num_orb(in2)
         write(unit,*) "iatom ",iatom," ineigh ",ineigh, &
                      " jatom ",jatom," in1 ",in1," in2 ",in2, &
                      " nmu ",nmu," nnu ",nnu," ioff ",ioff," joff ",joff 
         do inu = 1, nnu
          jnu = inu + joff
          do imu = 1, nmu
           jmu = imu + degelec(iatom)
           !write(*,'(A,e20.10,e20.10,A)',advance='no') M(imu,inu,ineigh,iatom)
           write(unit,'(e20.10)',advance='no') real(M(imu,inu,ineigh,iatom))
          end do ! do inu
          write(unit,*) 
         end do ! do imu
        end do ! do ineigh
    end do ! do iatom
    if (is_file_output) then
       close(unit)
    else
         write(*,*) "debug_writeBlockedMat(name= ", name, " ) DONE"
         write(*,*) "  "
    end if
end subroutine debug_writeBlockedMat

subroutine debug_writeIntegral( ifile, interaction, isub, in1, in2, index )
    use integrals
    !xintegral_2c(non2c,1,jxx,in1,in2)
    !xintegral_2c(integral,ipoint,itype,in1,in2)
    !      interpolate_1d (interaction, isub,  in1,  in2, non2c, ioption, xin,     yout,         dfdx         )
    !call interpolate_1d (interaction, ideriv, in1, in2, index, iforce, distance, slist(index), dslist(index))
    !ndex_coulomb = nssh(in1)*nssh(in2)
    !interaction = 12
    !ideriv = 0
    !do index = 1, index_coulomb
    ! call interpolate_1d (interaction, ideriv, in1, in2, index, iforce, distance, slist(index), dslist(index))
    !end do
    ! --- from read_2c()
    !itype = ind2c(interaction,isorp)
    !z2cmax(itype,in1,in2) = zmax
    !numz2c(itype,in1,in2) = numz
    ! ===== Parameters
    !character(len=*), intent(in) :: fname
    integer, intent(in) :: ifile
    integer, intent(in) :: interaction
    integer, intent(in) :: isub
    integer, intent(in) :: in1
    integer, intent(in) :: in2
    integer, intent(in) :: index
    ! ===== Variables
    integer itype, numz,i
    real    zmax,dz, val
    ! ===== Body
    itype = ind2c(interaction,isub)
    !nnum = numz2c(jxx,in1,in2)
    numz  = numz2c(itype,in1,in2)
    zmax  = z2cmax(itype,in1,in2)
    dz    = zmax/numz
    !open( 164646, file=fname, status='unknown' ) 
    !write(*,*) "debug_writeIntegral() shape(xintegral_2c) ",  shape(xintegral_2c)
    !write(*,*) "debug_writeIntegral() shape(splineint_2c) ",  shape(splineint_2c)
    !write(*,*) "debug_writeIntegral() index,itype,in1,in2, numz ",  index,itype,in1,in2,   numz
    !val = xintegral_2c(1,1,1,1,1)
    !write(*,*) "debug_writeIntegral() val ", val
    do i=1,numz
        !write(*,*) "debug_writeIntegral() i ", i
        !write(*,*) "debug_writeIntegral() index,i,itype,in1,in2 ", index,i,itype,in1,in2
        !splineint_2c
        write (ifile,*) dz*i, splineint_2c(1,index,i,itype,in1,in2)
        !write (ifile,*) dz*i, xintegral_2c(index,i,itype,in1,in2)
        !write (*,*) dz*i, xintegral_2c(index,i,itype,in1,in2)
        !write (*,*) xintegral_2c(index,i,itype,in1,in2)
    end do ! numz
    !close(164646)
end subroutine debug_writeIntegral


subroutine debug_writeIntegralSet( fname, interaction, isub, in1, in2 )
    use interactions
    ! from    doscentros (interaction, isub, iforce, in1, in2, in3,  distance, eps, deps, sx, spx) 
    !do index = 1, index_max2c(in1,in3)
    !    if ( switch ) then
    !        call interpolate_1d (interaction, isub, in1, in2, index, iforce, distance, slist(index), dslist(index))
    !    else
    !        call interpolate_1d (interaction, isub, in1, in3, index, iforce, distance, slist(index), dslist(index))
    !    end if
    !end do
    ! ===== Parameters
    character(len=*), intent(in) :: fname
    integer, intent(in) :: interaction
    integer, intent(in) :: isub
    integer, intent(in) :: in1
    integer, intent(in) :: in2
    ! ===== Variables
    integer ifile, index
    ! ===== Body
    ifile = 164646
    write (ifile,*) " debug_writeIntegralSet ", fname
    open( ifile, file=fname, status='unknown' ) 
    do index = 1, index_max2c(in1,in2)
        ! call interpolate_1d (interaction, isub, in1, in2, index, iforce, distance, slist(index), dslist(index))
        write (ifile,*) " ##### index ", index
        call debug_writeIntegral( ifile, interaction, isub, in1, in2, index )
    end do ! index
    close(ifile)
end subroutine debug_writeIntegralSet

end module debug