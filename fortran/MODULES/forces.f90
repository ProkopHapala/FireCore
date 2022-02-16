        module forces
!$ volatile dewald, fewald, flrew, dipp, sp_mat, spVNL, tp_mat, dusr, dxcv 
!$ volatile fro, ft, ftot, fana, fanl, faxc, fotna, fotnl, fotxc, f3naa 
!$ volatile f3nab, f3nac, f3nla, f3nlb, f3nlc, f3xca, f3xcb, f3xcc, faca 
!$ volatile faxc_ca, fotca, fotxc_ca, f3caa, f3cab, f3cac, f3xca_ca, f3xcb_ca
!$ volatile f3xcc_ca, fcoulomb, fxcnu

!CHROM
!		procedure (), pointer :: getforces
!END CHROM

! Ewald and long-range forces
         real, dimension (:, :, :), allocatable :: dewald
         real, dimension (:, :),    allocatable :: fewald
         real, dimension (:, :),    allocatable :: flrew

!        real, dimension (:, :), allocatable :: flrew_qmmm     ! IF_DEF_QMMM_END
!        real, dimension (:, :), allocatable :: ftot_dftd3     ! IF_DEF_DFTD3_END
 
! Derivatives of interactions
         real, dimension (:, :, :, :, :),    allocatable :: dipp
         real, dimension (:, :, :, :),       allocatable :: dippcm
         real, dimension (:, :, :, :, :, :), allocatable :: dippc
         real, dimension (:, :, :, :, :),    allocatable :: sp_mat
         real, dimension (:, :, :, :, :),    allocatable :: spm_mat
         real, dimension (:, :, :, :, :),    allocatable :: spVNL
         real, dimension (:, :, :, :, :),    allocatable :: tp_mat

! Forces
         real, dimension (:, :), allocatable, target :: ftot
         real, dimension (:, :), allocatable :: ftotold
         real, dimension (:, :), allocatable :: ftotnew

         real, dimension (:, :), allocatable :: dusr
         real, dimension (:, :), allocatable :: dxcv
         real, dimension (:, :), allocatable :: fro
         real, dimension (:, :), allocatable :: ft

         real, dimension (:, :), allocatable :: deltaF

         real, dimension (:, :, :), allocatable :: fana
         real, dimension (:, :, :), allocatable :: fanl
         real, dimension (:, :, :), allocatable :: faxc
         real, dimension (:, :, :), allocatable :: fotna
         real, dimension (:, :, :), allocatable :: fotnl
         real, dimension (:, :, :), allocatable :: fotxc

         real, dimension (:, :), allocatable :: f3naa, f3nab, f3nac
         real, dimension (:, :), allocatable :: f3nla, f3nlb, f3nlc
         real, dimension (:, :), allocatable :: f3xca, f3xcb, f3xcc

! Forces - OSLXC
         real, dimension (:, :, :), allocatable :: dxcdcc
         real, dimension (:, :, :, :, :), allocatable :: arhop_off
         real, dimension (:, :, :, :, :), allocatable :: arhopij_off 
         real, dimension (:, :, :, :, :), allocatable :: rhop_off
         real, dimension (:, :, :, :, :), allocatable :: rhopij_off
         real, dimension (:, :, :, :, :), allocatable :: rhop_on
         real, dimension (:, :, :, :, :), allocatable :: arhop_on

! Forces - DOGS
         real, dimension (:, :, :), allocatable :: faca
         real, dimension (:, :, :), allocatable :: faxc_ca
         real, dimension (:, :, :), allocatable :: fotca
         real, dimension (:, :, :), allocatable :: fotxc_ca

         real, dimension (:, :), allocatable :: f3caa, f3cab, f3cac
         real, dimension (:, :), allocatable :: f3xca_ca, f3xcb_ca, f3xcc_ca

!         real, dimension (:, :), allocatable :: fcoulomb     ! IF_DEF_Hubbard_END
!         real, dimension (:, :), allocatable :: fxcnu        ! IF_DEF_Hubbard_END

!         real, dimension (:, :, :), allocatable :: dxcdcc_zw     ! IF_DEF_ZW_END

!         real, dimension (:, :), allocatable :: fvdw          ! IF_DEF_vdW_END
!         real, dimension (:, :), allocatable :: fharmonic     ! IF_DEF_harmonic_END
!         real, dimension (:, :, :), allocatable :: fxcro      ! IF_DEF_Gauss_END
!         real, dimension (:, :), allocatable :: fbias         ! IF_DEF_bias_END

        end module
