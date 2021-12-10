
module debug
    implicit none
    contains

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

subroutine debug_writeBlockedMat( name, M )
    use configuration
    use neighbor_map
    use interactions
    ! allocate (s_mat (numorb_max, numorb_max, neigh_max, natoms))
    ! ===== Parameters
    character(len=*)                                               , intent(in)  :: name
    real, dimension( numorb_max, numorb_max, neigh_max, natoms ), intent (in) :: M
    ! ===== Variables
    integer iatom, jatom, ineigh, inu, imu,  jnu, jmu, in1,in2
    ! ===== Body
    open( 164646, file=name, status='unknown' ) 
    do iatom = 1, natoms
        in1 = imass(iatom)
        do ineigh = 1, neighn(iatom)
         !mbeta = neigh_b(ineigh,iatom)
         jatom = neigh_j(ineigh,iatom)
         in2 = imass(jatom)
         write(164646,*) "iatom ",iatom," ineigh ",ineigh 
         do inu = 1, num_orb(in2)
          jnu = inu + degelec(jatom)
          do imu = 1, num_orb(in1)
           jmu = imu + degelec(iatom)
           !write(*,'(A,e20.10,e20.10,A)',advance='no') M(imu,inu,ineigh,iatom)
           write(164646,'(e20.10)',advance='no') real(M(imu,inu,ineigh,iatom))
          end do ! do inu
          write(164646,*) 
         end do ! do imu
        end do ! do ineigh
    end do ! do iatom
    close(164646)
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