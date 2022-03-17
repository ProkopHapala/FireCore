        module density

!$ volatile bbnkre, bbnkim, blowim, blowre, cape, rho

! These arrays contain the coefficients to calculate the density matrix
! and the Lowdin charges:
         real, dimension (:, :, :), allocatable, target :: bbnkre
         real, dimension (:, :, :), allocatable, target :: bbnkim
         real, dimension (:, :, :), allocatable, target :: blowim
         real, dimension (:, :, :), allocatable, target :: blowre
         real, dimension (:, :), allocatable, target :: eigen_k

! These arrays store the densities.
         real, dimension (:, :, :, :), allocatable :: cape
         real, dimension (:, :, :, :), allocatable :: rho
         real, dimension (:, :, :, :), allocatable :: rhoPP

! These arrays store the densities in excited state.
!         real, dimension (:, :, :, :), allocatable :: cape_es
!         real, dimension (:, :, :, :), allocatable :: rho_es
!         real, dimension (:, :, :, :), allocatable :: rhoPP_es

! IF_DEF_GRID
          real, dimension (:, :), allocatable :: rhoA
          real, dimension (:, :, :, :), allocatable :: rho_old
! END_DEF_GRID

! These arrays store stuff related to the average density.
! Used in the OLSXC exchange-correlation interactions.
         real, dimension (:, :, :   ), allocatable :: rho_on
         real, dimension (:, :, :   ), allocatable :: rhoi_on
         real, dimension (:, :, :, :), allocatable :: rho_off
         real, dimension (:, :, :, :), allocatable :: rhoij_off
         real, dimension (:, :, :   ), allocatable :: arho_on
         real, dimension (:, :, :   ), allocatable :: arhoi_on
         real, dimension (:, :, :, :), allocatable :: arho_off
         real, dimension (:, :, :, :), allocatable :: arhoij_off
! JOM-nonadiabatic
        !  real, dimension ( :, :), allocatable :: foccupy_na
        !  integer, dimension ( :, :), allocatable :: ioccupy_na
! VLADA-nonadiabatic & cdft-mdet
        !  real, dimension (:, :, :), allocatable :: bbnkre_o 
        !  real, dimension (:, :, :), allocatable :: blowre_o  
        !  real, dimension (:, :, :), allocatable :: blowre_hist
        !  real, dimension ( :, :), allocatable :: foccupy_wr_gs    
        !  real, dimension ( :, :), allocatable :: foccupy_wr_gs_check  
        !  real, dimension ( :, :), allocatable :: foccupy_wr_proj
        !  real, dimension ( :, :), allocatable :: foccupy_na_o
        !  integer, dimension ( :, :), allocatable :: ioccupy_na_o
       end module density
