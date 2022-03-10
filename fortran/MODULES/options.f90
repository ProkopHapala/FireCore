 module options
! This module provide information about basic options of the Fireball task
! ===========================================================================

! ---------------------------------------------------------------------------
! Toggles read from the options.input file.
! ---------------------------------------------------------------------------


 !       integer nprocs   iordern  ! IF_DEF_ORDERN_END

        integer iparam_file

! ====== Electronic Theory Modes

        ! Level of theory (0=> harris, 1=> idogs, 2=> extended hubbard)
        integer itheory
! Level of exchange-correlation theory
! 0 => Horsfield,
! 1 => Generalized Sankey-Niklewski (GSN)
! 2 => McWEDA
        integer itheory_xc

        integer igsn         ! Generalized Snakye Niklewsky (even simpler than Harris?)
        integer iharris      ! simplest non-SCF ( SNXC Snakye-Niklewski )
        integer idogs        ! SCF on top of Harris
        integer imcweda      ! default SCF theory, OSLXC
!        integer ixczw                  ! Second order XC theory ( Diego )                               ! IF_DEF_IXCZW_END
!        integer ihorsfield  ! version of SCF SN not work well for condensed systems, McWeda is corrected version of this  ! IF_DEF_HORSFIELD_END
!        integer ihubbard    ! simple kind of Hubbard_U / LDA+U, not used, can remove                    ! IF_DEF_HUBBARD_END
        integer iks          ! Kohn-Sahm Grid                                                            ! IF_DEF_GRID_END
        integer igrid        ! the grid projection                                                       ! IF_DEF_GRID_END
        integer iwrtewf      ! write out wavefunctions                                                   ! IF_DEF_GRID_END
!        integer igauss     ! Use gaussian approximation to 3-center                                     ! IF_DEF_GAUSS_END 

        integer idipole     ! Long range term with XYZ dipole (correction of long range electrostatics, problem in PBC)
!       integer V_intra_dip            ! Intra-atomic dipolar potential   ( Diego )                                     ! IF_DEF_IXCZW_END

! ======  Optional Modules - Electronic Theory

!        integer iephc                  ! electron-phonon coupling                          ! IF_DEF_ephc_END
!        integer itdse                  ! Time Propagatin of wavefunction                   ! IF_DEF_TDSE_END
!        integer imdet                  ! Tully Minimal Switches algortithm  JOM-add        ! IF_DEF_MDET_END
!        integer iProjWF                ! do projection within MDET simulations             ! IF_DEF_MDET_END

!        integer ipathintegral          ! Add quantum effects - path integral               ! IF_DEF_pathintegral_END
!        integer iordern                ! Perform linear-scaling algorithm                  ! IF_DEF_ordern_END
!        integer itrans                 ! do transport calculations                         ! IF_DEF_trans_END
!        integer ibias                  ! to apply the potential ramp (bias voltage) on hamiltonian ! IF_DEF_bias_END

!        integer igap                   ! to introduce LDA gap corrections GAP ENRIQUE-FF   ! IF_DEF_gap_END
!        integer icDFT                  ! do the constrain DFT                              ! IF_DEF_CDFT_END

!        integer iqmmm                  ! QM/MM Electrostatic embedding                         ! IF_DEF_QMMM_END
!        integer mix_embedding          ! mix electrostatic and mechanical embedding in QM/MM   ! IF_DEF_QMMM_END
!        real    cut_embedding          ! cutoff for the mix embedding                          ! IF_DEF_QMMM_END

! ====== Optional Modules - Ionic / Classical

!        integer idynmat                ! Dynamical matrix simulation               ! IF_DEF_dynmat_END
!        integer iharmonic              ! whether to attach harmonic oscillators    ! IF_DEF_harmonic_END

!        integer iumbrella              ! Do umbrella sampling (on distances)       ! IF_DEF_umberla_END
!        integer ithermoint             ! do thermodynamic integration              ! IF_DEF_thermoint_END
!        integer ireducekpts            ! whether to reduce kpts by atm symmetry    ! IF_DEF_autoK_END   ! IF_DEF_PBC_END
!        integer ineb                   ! do Nudged Elastic Band Method             ! IF_DEF_NEB_END
!        integer iclassicMD             ! classical forcefields  (Zdenka Chromcova) ! IF_DEF_classicMD_END

!        integer ivdw                   ! Include van der Waals interactions        ! IF_DEF_VDW_END
!        integer idftd3                 ! DFTD3 corrections                         ! IF_DEF_DFTD3_END
! ! IF_DEF_DFTD3
!         real  dftd3_s6, dftd3_rs6, dftd3_s18, dftd3_rs18, dftd3_alp
!         character (len = 40) dftd3_func
!         integer dftd3_version
!         logical dftd3_tz 
!         real dftd3_params (5)
! ! END_DEF_DFTD3

! ====== Optional Modules - I/O & Misc

!        integer isocket                ! socket for i-pi ( Jesus )      ! IF_DEF_socket_END

! ====== Switches affecting Electronic Solution 

        integer ivec_2c   !   use vectorized versersion of interpolation rutines for 2-center integrals ?
        integer ivec_3c   !   use vectorized versersion of interpolation rutines for 3-center integrals ?
        integer i2dlin    !   use be-linear (rather than bi-cubic) interpolation for 3-center integrals ? 

!        integer imix                   ! which density mixer to use
        integer iforce                 ! Calculate forces ?                     ! IF_DEF_force_END
        integer icluster               ! Calculate gas-phase                    ! IF_DEF_PBC_END
        integer iimage                 ! How often we image the box.            ! IF_DEF_PBC_END

        integer iqout                  ! Charges to use (Lowdin or Mulliken)
!        integer ispin                  ! Spin interaction                      ! IF_DEF_spin_END

        integer ifixcharge             ! Fix the charges                        ! IF_DEF_SCF_END
        integer ifixneigh              ! Fix the neighbors                      ! IF_DEF_neighbors_END

! ====== Switches affecting Ionic Solution / Classical

!        integer iquench                ! Quenching option                       ! IF_DEF_move_END
!        integer iensemble              ! Which ensemble?                        ! IF_DEF_move_END
!        integer iendtemp               ! toggles T_final for MD calculations    ! IF_DEF_move_END

! ======= Debug / Test / Log

        integer timing_verbosity
        integer idebugWrite
        integer verbosity       ! 0 minimal, 10 max         
        integer ntpr            ! Every ntpr steps, information will be printed , if = 0 nothing will be write until the end        
        integer restartxyz      ! if = 1, start the simulation from restart.xyz
        integer inputxyz        ! if = 1, the coordinate input file is a .xyz

        integer itestrange
        integer iconstraints (4)  ! big (4) constraints
        integer ioff2c (1:24)    ! for diagnostic purposes
        integer ioff3c (1:4)
        real    testrange

        integer iwrtxyz



 end module options
