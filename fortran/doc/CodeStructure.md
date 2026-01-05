
## Electronic Problem Options







* Core ?
    * `itheory`
    * `itheory_xc`
* Optional
    * `iks`         Kohn-Sham grid
    * `iqmmm`       Quantum-Mechanics/Molecular-Mechanins
    * `iordern`     Linear-Scaling solver  (O(n))
    * `ibias`       Extarnal potential ramp (e.g. Bias Voltage)
    * `imdet`       MDET 
    * iCDFT
* Possibly remove ?
    * `iforce`      Derivatives of Energy(Hamiltonian,Den-Mat,S-Mat) to obtain forces on atoms
* Not sure 
    * `idipole`    Dipole correction; ( seems to be connected with QMMM )
    * `igauss`     Calculate matrix elements (H-mat,S-mat) using Gaussian functions
    * `ivdw`       van-der-Waals correction
    * `iharmonic`  Harmonic correction to energy ? ( Zero energy vibrations ? ) 
* Minor / IO
    * `iwrtfpieces`  Write individual components of Hamiltonian (?) (For debugging?)

## Notes
 * `(D)assemble_ca_*`        ? DOGS ?  What is DOGS ? 
 * `(D)assemble_eh_*`        (?) Extended Hubbard
 * `(D)assemble_lr_*`        long-range ewald 
 * `(D)assemble_hxc_*`       Horsfield exchange-correlation
 * `(D)assemble_snsxc_*`     Sankey-Niklewski two-center  exchange correlation
   * `(D)assemble_olsxc_*`   What is the difference from  `snxc`  ? 
 * `(D)assemble_scissor`        Koopman correction potential 


## NOTES:
 * `/ASSEMBLERS/getenergy.f90` contains nice overview of different electronics structure options (`itheory` etc.)


## COMMENTING MODULES

* `END_DEF_*`
* `IF_DEF_*`
* `IF_DEF_*_END`

* `*_hxc`  is Horsfield

## All options:

#### Electronic Theory Modes

* `itheory`    Level of theory (0=> harris; 1=> idogs; 2=> extended hubbard, 3=> Kohn-Sham )

* `itheory_xc`  Level of exchange-correlation theory ( 0=> Horsfield; 1=>Generalized Sankey-Niklewski (GSN); 2=> McWEDA, 4=> ZW)
* `igsn`        Generalized Snakye Niklewsky (even simpler than Harris?)
* `iharris`     implest non-SCF ( SNXC Snakye-Niklewski )
* `idogs`       SCF on top of Harris
* `imcweda`     default SCF theory, OSLXC
* `ixczw`       Second order XC theory ( Diego )                               ! IF_DEF_IXCZW_END
* `ihorsfield`  version of SCF SN not work well for condensed systems, McWeda is corrected version of this  ! IF_DEF_HORSFIELD_END
* `ihubbard`    simple kind of Hubbard_U / LDA+U, not used, can remove                    ! IF_DEF_HUBBARD_END
* `iks`         Kohn-Sahm Grid                                                            ! IF_DEF_KS_END

* `igauss`      Use gaussian approximation to 3-center                                    ! IF_DEF_GAUSS_END 

* `idipole`     Long range term with XYZ dipole (correction of long range electrostatics, problem in PBC)
* `V_intra_dip` Intra-atomic dipolar potential   ( Diego )                                     ! IF_DEF_IXCZW_END

####  Optional Modules - Electronic Theory

* `iephc`            electron-phonon coupling                          ! IF_DEF_ephc_END
* `itdse`            Time Propagatin of wavefunction                   ! IF_DEF_TDSE_END
* `imdet`            Tully Minimal Switches algortithm  JOM-add        ! IF_DEF_MDET_END
* `iProjWF`          do projection within MDET simulations             ! IF_DEF_MDET_END

* `ipathintegral`    Add quantum effects - path integral               ! IF_DEF_pathintegral_END
* `iordern`          Perform linear-scaling algorithm                  ! IF_DEF_ordern_END
* `itrans`           do transport calculations                         ! IF_DEF_trans_END
* `ibias`            to apply the potential ramp (bias voltage) on hamiltonian ! IF_DEF_bias_END

* `igap`             to introduce LDA gap corrections GAP ENRIQUE-FF   ! IF_DEF_gap_END
* `icDFT`            do the constrain DFT                              ! IF_DEF_CDFT_END

* `iqmmm`            QM/MM Electrostatic embedding                         ! IF_DEF_QMMM_END
* `mix_embedding`    mix electrostatic and mechanical embedding in QM/MM   ! IF_DEF_QMMM_END

#### Optional Modules - Ionic / Classical

* `idynmat`      Dynamical matrix simulation               ! IF_DEF_dynmat_END
* `iharmonic`    whether to attach harmonic oscillators    ! IF_DEF_harmonic_END

* `iumbrella`    Do umbrella sampling (on distances)       ! IF_DEF_umberla_END
* `ithermoint`   do thermodynamic integration              ! IF_DEF_thermoint_END
* `ireducekpts`  whether to reduce kpts by atm symmetry    ! IF_DEF_autoK_END   ! IF_DEF_PBC_END
* `ineb`         do Nudged Elastic Band Method             ! IF_DEF_NEB_END
* `iclassicMD`   classical forcefields  (Zdenka Chromcova) ! IF_DEF_classicMD_END

* `ivdw`         Include van der Waals interactions        ! IF_DEF_VDW_END
* `idftd3`       DFTD3 corrections                         ! IF_DEF_DFTD3_END


#### Optional Modules - I/O & Misc

* `igrid`    the grid projection         ! IF_DEF_KS_END
* `isocket`  socket for i-pi ( Jesus )   ! IF_DEF_socket_END

#### Switches affecting Electronic Solution 

* `imix`     which density mixer to use
* `iforce`   Calculate forces ?                  ! IF_DEF_force_END
* `icluster` Calculate gas-phase                 ! IF_DEF_PBC_END
* `iimage`   How often we image the box.         ! IF_DEF_PBC_END

* `iqout`    Charges to use (Lowdin or Mulliken)
* `ispin`    Spin interaction                    ! IF_DEF_spin_END

* `ifixcharge`  Fix the charges                  ! IF_DEF_SCF_END
* `ifixneigh`   Fix the neighbors                ! IF_DEF_neighbors_END

#### Switches affecting Ionic Solution / Classical

* `iquench`    Quenching option                     ! IF_DEF_move_END
* `iensemble`  Which ensemble?                      ! IF_DEF_move_END
* `iendtemp`   toggles T_final for MD calculations  ! IF_DEF_move_END

#### Debug / Test / Log

* `verbosity`   0 minimal, 10 max         
* `ntpr`        Every ntpr steps, information will be printed , if = 0 nothing will be write until the end        
* `restartxyz`  if = 1, start the simulation from restart.xyz
* `inputxyz`    if = 1, the coordinate input file is a .xyz





