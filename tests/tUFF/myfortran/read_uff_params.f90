subroutine read_uff_params
  use uff
  implicit none
  real(8), parameter :: pi = acos(-1.d0)
  real(8), parameter :: au_to_Ang = 0.52917d0
  real(8), parameter :: eV_to_kcal = 1.602176634d-19 / 4184.d0 * 6.02214076d23
  integer :: i
  
  open(unit=10,file='params.txt')
  read(10,*)
  read(10,*) nufftypes
  allocate(sym(nufftypes))
  allocate(r_I(nufftypes))     ! Ang
  allocate(theta_0(nufftypes)) ! deg
  allocate(x_I(nufftypes))     ! Ang
  allocate(D_I(nufftypes))     ! kcal/mol
  allocate(zeta(nufftypes))    ! .
  allocate(Zstar_I(nufftypes)) ! e
  allocate(chi(nufftypes))     ! eV
  allocate(J_I(nufftypes))       ! eV
  allocate(R(nufftypes))       ! Ang
  allocate(eta(nufftypes))     ! 1/au
  allocate(V_I(nufftypes))     ! kcal/mol
  allocate(U_I(nufftypes))     ! kcal/mol
  allocate(chi_reax(nufftypes))   ! eV
  allocate(eta_reax(nufftypes))   ! eV
  allocate(gamma_reax(nufftypes)) ! 1/Ang
  do i = 1, nufftypes
     read(10,*) sym(i), r_I(i), theta_0(i), x_I(i), D_I(i), zeta(i), Zstar_I(i), chi(i), J_I(i), R(i), eta(i), V_I(i), U_I(i), chi_reax(i), eta_reax(i), gamma_reax(i)
     theta_0(i) = theta_0(i) / 180.d0 * pi ! rad
     chi(i) = chi(i) * eV_to_kcal ! kcal/mol
     J_I(i) = J_I(i) * eV_to_kcal ! kcal/mol
     eta(i) = eta(i) / au_to_Ang ! 1/Ang
  end do
  close(10)

  return
end subroutine read_uff_params
