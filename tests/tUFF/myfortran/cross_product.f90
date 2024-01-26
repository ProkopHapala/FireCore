subroutine cross_product ( a, b, n )
  implicit none
  real(8), intent(in) :: a(3), b(3)
  real(8), intent(out) :: n(3)

  n(1) = a(2) * b(3) - a(3) * b(2)
  n(2) = a(3) * b(1) - a(1) * b(3)
  n(3) = a(1) * b(2) - a(2) * b(1)

  return
end subroutine cross_product
