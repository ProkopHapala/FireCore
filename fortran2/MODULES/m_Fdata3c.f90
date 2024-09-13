module m_Fdata3c

    implicit none

    type, public :: Fdata3c
        integer nx, ny, nME
        real    xmax, ymax
        real    hx, hy
        !allocate (bcna_vec (5,numXmax, numYmax, ME3c_max, 0:isorpmax, nspecies**3))
        real,allocatable :: data(:,:,:,:)
        contains
        !procedure  :: new  => t_data3c_new
        !procedure  :: set  => t_data3c_set
        !procedure  :: load => t_data3c_load
        !procedure  :: interpolate_2d => t_data3c_interpolate_2d
    end type

    !public :: tt_data3c_load !(this,xin, yin, iforce, Q_L, dQ_Ldx, dQ_Ldy )

    contains

    subroutine fdata3c_new(this, nx,ny,nME)
        implicit none
        type    ( Fdata3c ) :: this
        integer,intent(in)   :: nx, ny, nME
        this%nx=nx
        this%nx=nx
        this%nME=nME
        allocate( this%data(5,nx,ny,nME) )  
    end subroutine

    subroutine fdata3c_set(this, hx, hy, xmax, ymax )
        implicit none
        type( Fdata3c ) :: this
        real,intent(in)  :: hx, hy, xmax, ymax
        this%hx=hx
        this%hy=hy
        this%xmax=xmax 
        this%ymax=ymax
    end subroutine

    ! subroutine t_data3c_load(this, isorp,index, d1,d2,d3,d4,d5 )
    !     use dimensions
    !     !use integrals
    !     use interactions
    !     use configuration, only : nspecies
    !     implicit none
    !     type   ( t_data3c ),intent(inout) :: this
    !     integer,intent(in)  :: isorp,index
    !     real   ,intent(in), dimension( numXmax, numYmax, ME3c_max, 0:isorpmax, nspecies**3)  :: d1,d2,d3,d4,d5
    !     integer i
    !     do i=1,this%nME
    !         this%data(1,i,:,:) = d1( :,:, i, isorp,index )
    !         this%data(2,i,:,:) = d2( :,:, i, isorp,index )
    !         this%data(3,i,:,:) = d3( :,:, i, isorp,index )
    !         this%data(4,i,:,:) = d4( :,:, i, isorp,index )
    !         this%data(5,i,:,:) = d5( :,:, i, isorp,index )
    !     end do
    ! end subroutine

    ! subroutine t_data3c_interpolate_2d(this,xin, yin, iforce, Q_L, dQ_Ldx, dQ_Ldy )
    !     use dimensions
    !     use interactions
    !     use constants_fireball
    !     use timing
    !     implicit none
    ! ! ============= Arguments
    !     type( t_data3c )  :: this
    !     real, intent (in) :: xin
    !     real, intent (in) :: yin
    !     integer, intent (in) :: iforce
    !     real, intent (out), dimension (5,numorb_max) ::  Q_L   
    !     real, intent (out), dimension (5,numorb_max) :: dQ_Ldx 
    !     real, intent (out), dimension (5,numorb_max) :: dQ_Ldy
    ! ! ====== VARIABLES
    !     integer imidx, imidy
    !     integer k, ik
    !     real px, py
    !     real,dimension (5,this%nME)   :: f1m1, f0p3, f1p3, f2p1
    !     real,dimension (5,this%nME)   :: bb0,bb1,bb2,bb3
    !     real,dimension (5,this%nME,4) :: g,gp
    !     real inv_hx
    !     real inv_hy
    ! ! =========== BODY
    !     !ncall_interpolate_2d_vec=ncall_interpolate_2d_vec+1

    !     inv_hx = 1/this%hx
    !     inv_hy = 1/this%hy

    !     imidx = int(xin*inv_hx) + 1
    !     if (imidx .lt. 2) then
    !       imidx = 2
    !     else if (imidx .gt. this%nx - 2) then
    !       imidx = this%nx - 2
    !     end if

    !     imidy = int(yin*inv_hy) + 1
    !     if (imidy .lt. 2) then
    !       imidy = 2
    !     else if (imidy .gt. this%ny - 2) then
    !       imidy = this%ny - 2
    !     end if

    !     px = xin*inv_hx - (imidx - 1)
    !     py = yin*inv_hx - (imidy - 1)
 
    !     do k = 1, 4
    !         ik = imidx + k-2
    !         f1m1(:,:) =   this%data(:,:,ik,imidy - 1 )
    !         f0p3(:,:) = 3*this%data(:,:,ik,imidy     )
    !         f1p3(:,:) = 3*this%data(:,:,ik,imidy + 1 )
    !         f2p1(:,:) =   this%data(:,:,ik,imidy + 2 )
    !         bb3 = -   f1m1 +   f0p3 -   f1p3 + f2p1
    !         bb2 =   3*f1m1 - 2*f0p3 +   f1p3
    !         bb1 = - 2*f1m1 -   f0p3 + 2*f1p3 - f2p1
    !         bb0 =   2*f0p3
    !         g(:,:,k)                     = ((  bb3*py +   bb2)*py + bb1)*py + bb0
    !         if (iforce .eq. 1) gp(:,:,k) = ((3*bb3*py + 2*bb2)*py + bb1)
    !     end do ! k
    !     f1m1(:,:) =   g(:,:,1)
    !     f0p3(:,:) = 3*g(:,:,2)
    !     f1p3(:,:) = 3*g(:,:,3)
    !     f2p1(:,:) =   g(:,:,4) 
    !     bb3 = -   f1m1 +   f0p3 -   f1p3 + f2p1
    !     bb2 =   3*f1m1 - 2*f0p3 +   f1p3
    !     bb1 = - 2*f1m1 -   f0p3 + 2*f1p3 - f2p1
    !     bb0 =   2*f0p3
    !     Q_L(:,:this%nME) = (((bb3*px + bb2)*px + bb1)*px + bb0)/36
    !     if (iforce .eq. 1) then
    !         dQ_Ldx(:,:this%nME) = ((3*bb3*px + 2*bb2)*px + bb1)*inv_hx/36
    !         f1m1(:,:) =   gp(:,:,1)
    !         f0p3(:,:) = 3*gp(:,:,2)
    !         f1p3(:,:) = 3*gp(:,:,3)
    !         f2p1(:,:) =   gp(:,:,4) 
    !         bb3 = -   f1m1 +   f0p3 -   f1p3 + f2p1
    !         bb2 =   3*f1m1 - 2*f0p3 +   f1p3
    !         bb1 = - 2*f1m1 -   f0p3 + 2*f1p3 - f2p1
    !         bb0 =   2*f0p3
    !         dQ_Ldy(:,:this%nME) = (((bb3*px + bb2)*px + bb1)*px + bb0)*inv_hy/36
    !     end if
    ! return
    ! end

end module Fdata3c