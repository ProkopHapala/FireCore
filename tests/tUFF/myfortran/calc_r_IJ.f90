subroutine calc_r_IJ ( i, j, BO, rIJ )
  use uff, only: r_I, chi
  implicit none
  integer, intent(in) :: i, j
  real(8), intent(in) :: BO
  real(8), intent(out) :: rIJ
  real(8) :: rBO, rEN

  rBO = -0.1332d0 * ( r_I(i) + r_I(j) ) * log( BO )
  
  rEN = r_I(i) * r_I(j) * ( sqrt(chi(i)) - sqrt(chi(j)) )**2 / ( chi(i) * r_I(i) + chi(j) * r_I(j) )
  
  rIJ = r_I(i) + r_I(j) + rBO - rEN

  return
end subroutine calc_r_IJ
