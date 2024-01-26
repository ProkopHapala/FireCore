subroutine calc_dist ( pos_i, pos_j, dist )
  use conf, only: side
  implicit none
  real(8), intent(in) :: pos_i(3), pos_j(3)
  real(8), intent(out) :: dist
  integer :: ix, iy, iz, k
  real(8) :: rawdist_frac(3), dist_frac(3), dist_cart(3), dist_current

  ! compute raw distance in fractional coordinates
  rawdist_frac(:) = pos_i(:) - pos_j(:)

  ! apply standard minimum image convention
  rawdist_frac(:) = rawdist_frac(:) - dble(nint(rawdist_frac(:)))

  ! loop over the 27 images to find the shortest Cartesian distance
  dist = 1.d99
  do ix = -1, 1
     dist_frac(1) = rawdist_frac(1) + dble(ix)
     do iy = -1, 1
        dist_frac(2) = rawdist_frac(2) + dble(iy)
        do iz = -1, 1
           dist_frac(3) = rawdist_frac(3) + dble(iz)
           do k = 1, 3
              dist_cart(k) = side(1,k)*dist_frac(1) + side(2,k)*dist_frac(2) + side(3,k)*dist_frac(3)
           end do
           dist_current = sqrt(dot_product(dist_cart,dist_cart))
           if ( dist_current < dist ) dist = dist_current
        end do
     end do
  end do

  return
end subroutine calc_dist
