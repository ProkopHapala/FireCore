subroutine calc_dist_eq ( i, j, dist_eq )
  use constants, only: radius_H, radius_C,  radius_N,  radius_O
  use conf, only: tipo
  implicit none
  integer, intent(in) :: i, j
  real(8), intent(out) :: dist_eq

  if ( tipo(i) == 'H' ) then
     dist_eq = radius_H
  else if ( tipo(i) == 'C' ) then
     dist_eq = radius_C
  else if ( tipo(i) == 'N' ) then
     dist_eq = radius_N
  else if ( tipo(i) == 'O' ) then
     dist_eq = radius_O
  end if

  if ( tipo(j) == 'H' ) then
     dist_eq = dist_eq + radius_H
  else if ( tipo(j) == 'C' ) then
     dist_eq = dist_eq + radius_C
  else if ( tipo(j) == 'N' ) then
     dist_eq = dist_eq + radius_N
  else if ( tipo(j) == 'O' ) then
     dist_eq = dist_eq + radius_O
  end if

  dist_eq = 0.5d0 * dist_eq
  
  return
end subroutine calc_dist_eq
