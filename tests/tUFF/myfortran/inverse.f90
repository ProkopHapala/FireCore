subroutine inverse ( aa, c )
  implicit none 
  integer, parameter :: n = 3
  real(8), intent(in) :: aa(n,n)
  real(8), intent(out) :: c(n,n)
  integer :: i, j, k
  real(8) :: a(n,n), L(n,n), U(n,n), b(n), d(n), x(n)
  real(8) :: coeff
  
  a(:,:) = aa(:,:)
  L(:,:) = 0.d0
  U(:,:) = 0.d0
  b(:) = 0.d0

  do k=1, n-1
     do i=k+1,n
        coeff=a(i,k)/a(k,k)
        L(i,k) = coeff
        do j=k+1,n
           a(i,j) = a(i,j)-coeff*a(k,j)
        end do
     end do
  end do

  do i=1,n
     L(i,i) = 1.d0
  end do
  do j=1,n
     do i=1,j
        U(i,j) = a(i,j)
     end do
  end do

  do k=1,n
     b(k)=1.d0
     d(1) = b(1)
     do i=2,n
        d(i)=b(i)
        do j=1,i-1
           d(i) = d(i) - L(i,j)*d(j)
        end do
     end do
     x(n)=d(n)/U(n,n)
     do i = n-1,1,-1
        x(i) = d(i)
        do j=n,i+1,-1
           x(i)=x(i)-U(i,j)*x(j)
        end do
        x(i) = x(i)/u(i,i)
     end do
     
     do i=1,n
        c(i,k) = x(i)
     end do
     b(k)=0.d0
  end do

  return
end subroutine inverse
