 module energy 
! This module defines variables of energies   
! ===========================================================================

! Energy and force contributions and terms 
        real atomic_energy
        real etot
        real ebs
        real etotnew
        real etotold
        real etotper
        real etotxc_1c
        real getot
        real getot_initial
        real getotper
        real deltaE
        real deltaFmax
        real uiiuee
        real uxcdcc
        real uxcdcc_sn
!        real uxcdcc_hf             ! IF_DEF_Horsfield_END
        real uxcdcc_ols
!        real uxcdcc_zw                       ! IF_DEF_ZW_END
!        real dc_v_intra_dip_1c               ! IF_DEF_ZW_END
!        real, dimension (3) :: duxcdcc_zw    ! IF_DEF_ZW_END
!        real uxcdcc_ks   ! IF_DEF_KS_END
!        real uhdcc_ks    ! IF_DEF_KS_END
! Ext Hubbard
!        real ehxcc       ! IF_DEF_Hubbard_END
!        real ehcoolc     ! IF_DEF_Hubbard_END
!        real Umuxc_1c     ! IF_DEF_Hubbard_END
!        real Uexc_1c      ! IF_DEF_Hubbard_END
!        real eqmmm        ! IF_DEF_QMMM_END
!        real etot_dftd3   ! IF_DEF_DFTD3_END

 end module energy 
