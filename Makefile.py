FORM_RHO = {
  '' : ['denmat','denmat_es','denmat_KS','fermie','denmata_ordern_fk','denmatb_ordern_fk','denmatc_ordern_fk','ss12_fk','build_rho',
        'build_rho_KS','build_Ji','koopman','build_nij','project_eh','project_wfmdet']
  'ORDERN' : ['chebft','denmata_ordern','denmatb_ordern','denmatc_ordern','denmat_fk','formrho_sparse','ss12']
}

KSPACE = {
  '' : ['kspace','kspace_KS','kspace_ordern_fk','diag_error','diag_k','diag_k_KS']
  'SCALAPACK' : ['blacsaba','kspace_MPI','kspace_MPI_slave','pclagetter','pclaputter']
  'ORDERN' : ['eandg','xeandg','formc_compact','formsh_compact','initguess','kspace_fk','kspace_ordern','kspace_ordern_init','kspace_ordern_slave',
        'ordern_init','qralg','set_dimensions','set_maxdimension','getsendrecv','getstepsize']
  'LAPACK95' : ['kspace_l95']
  'DOUBLE' : ['kspace2']
  'GAMMA' : ['kspace_withnok']
  'MPI-k' : ['kspace_MPI-k','diag_k-MPI']
}

INITMPI = {
  '' : ['mpi_declarations','init_MPI']
  'SCALAPACK' : []
  'LAPACK95' : []
}

UTIL = {
  '' : ['anderson','fixfrags','fixfrags2','hampiece','push_atoms','writeout_ac','writeout_cd','writeout_charges','writeout_dipole',
        'writeout_comph','writeout_neighbors','writeout_xv','writeout_neighborsPP','hamtrans','mixer','den2mesh','ew2mesh','ew2mesh_gamma',
        'postscf','den2mesh_import','broyden','louie','pulay','compute_neutral_rho_grid','ew2mesh_fourier','ew2mesh_kscan','ew2mesh_ARPES',
        'project_orb','project_orb_complex','excitations','kvaziband','band_project','ARPES_LCAO','writeout_eigenvec','writeout_eigenvec_G',
        'writeCoefsLCAO']
  'OPENMP' : ['anderson2']
  'DOUBLE' : ['anderson_l95']
  'GAMMA' : ['anderson2']
  'MPI-k' : ['anderson2','diag_k-MPI_slave','bcast_k-MPI','bcast_k-MPI_slave']
}

DASSEMBLERS = {
  '' : ['Dassemble_2c','Dassemble_3c','Dassemble_ca_2c','Dassemble_ca_3c','Dassemble_eh_2c','Dassemble_hxc_2c','Dassemble_hxc_3c',
        'Dassemble_lr','Dassemble_snxc_on','Dassemble_olsxc_on','Dassemble_olsxc_2c','Dassemble_olsxc_3c','Dassemble_snxc_2c','Dassemble_snxc_3c',
        'Dassemble_2c_PP','Dassemble_3c_PP','Dassemble_ca_olsxc_on','Dassemble_ca_snxc_on','Dassemble_ca_snxc_3c','Dassemble_ca_olsxc_3c',
        'Dassemble_ca_snxc_2c','Dassemble_ca_olsxc_2c','Dassemble_qmmm','Dassemble_qmmm_dip','Dassemble_ca_2c_dip','Dassemble_ca_3c_dip',
        'Dassemble_lr_dip','getforces_mcweda','getforces_eh','getforces_hxc','getforces_KS','getforces_classic','getforces_classic_RGL',
        'getforces_classic_vdw','getforces','getforces_classic_tersoff','getforces_zw','Dassemble_zw_2c_ct','Dassemble_zw_3c_ct','Dassemble_zw_2c_na',
        'Dassemble_zw_3c_na','Dassemble_zw_on_na']
  'OPENMP' : ['Dassemble_lr_OMP']
}

INTERACTIONS = {
  'SCALAPACK' : ['cl_value','Dtrescentros','doscentros','doscentrosPP','doscentrosDipY','doscentrosDipX','get_ewald','get_vdw','getHarmonic',
        'trescentros','unocentros','smoother','doscentrosS','trescentrosS','DtrescentrosS','dosgaussians','gelements_VXC','Dgelements_VXC',
        'gelementsG_VXC','DgelementsG_VXC','gelementsG_VNA','DgelementsG_VNA','gelementsGS_VXC','DgelementsGS_VXC','gelementsG_VNA_SH','DgelementsG_VNA_SH',
        'gelementsG_VNA_SH_RNA','DgelementsG_VNA_SH_RNA','trescentrosGHXC_VXC','DtrescentrosGHXC_VXC','trescentrosG_VNA','DtrescentrosG_VNA',
        'trescentrosG_VXC','DtrescentrosG_VXC','trescentrosGS_VXC','DtrescentrosGS_VXC','trescentrosG_VNA_SH','DtrescentrosG_VNA_SH','doscentrosGS_overlap',
        'gelementsGS_overlap','DgelementsGS_overlap','doscentrosG_overlap','gelementsG_overlap','DgelementsG_overlap','internalLambda',
        'tester2c']
  'ORDERN' : ['cl_value','Dtrescentros','doscentros','doscentrosPP','get_ewald_OMP','get_vdw','getHarmonic','trescentros','unocentros',
        'smoother','doscentrosS','trescentrosS','DtrescentrosS','dosgaussians','gelements_VXC','Dgelements_VXC','gelementsG_VXC','DgelementsG_VXC',
        'gelementsG_VNA','DgelementsG_VNA','gelementsGS_VXC','DgelementsGS_VXC','gelementsG_VNA_SH','DgelementsG_VNA_SH','gelementsG_VNA_SH_RNA',
        'DgelementsG_VNA_SH_RNA','trescentrosGHXC_VXC','DtrescentrosGHXC_VXC','trescentrosG_VNA','DtrescentrosG_VNA','trescentrosG_VXC',
        'DtrescentrosG_VXC','trescentrosGS_VXC','DtrescentrosGS_VXC','trescentrosG_VNA_SH','DtrescentrosG_VNA_SH','doscentrosGS_overlap',
        'gelementsGS_overlap','DgelementsGS_overlap','doscentrosG_overlap','gelementsG_overlap','DgelementsG_overlap','internalLambda',
        'tester2c']
}

ORDERN = {
  'ORDERN' : ['ordern','allocate_ordern','assemble_2c_ordern_final','assemble_2c_ordern_init','assemble_3c_ordern_final','assemble_ca_2c_ordern_final',
        'assemble_ca_3c_ordern_final','assemble_eh_2c_ordern_final','Dassemble_2c_ordern_final','Dassemble_3c_ordern_final','Dassemble_ca_2c_ordern_final',
        'Dassemble_ca_3c_ordern_final']
}

UTIL_SPARSE = {
  'ORDERN' : ['build_transpose','lanc','sparse_add','sparse_copy','sparse_getdimension','sparse_getpacksize','sparse_mask','sparse_mult',
        'sparse_norm2','sparse_pack','sparse_pack_elements','sparse_pack_indices','sparse_unpack','sparse_unpack_elements','sparse_unpack_indices',
        'sparse_vecprod']
}

ALLOCATIONS = {
  '' : ['allocate_f','allocate_h','allocate_neigh','allocate_rho','allocate_umb','allocate_steered','reallocate_f','reallocate_h',
        'reallocate_neigh','reallocate_rho','allocate_dos','allocate_grid','allocate_trans']
}

ASSEMBLERS = {
  '' : ['assemble_olsxc_1c','assemble_hxc_1c','assemble_2c','assemble_3c','assemble_ca_2c','assemble_3c_PP','assemble_ca_3c','assemble_eh_2c',
        'assemble_eh_usr','assemble_F','assemble_hxc_2c','assemble_hxc_3c','assemble_lr','assemble_sVNL','assemble_usr','buildh','assemble_olsxc_on',
        'assemble_olsxc_off','build_olsxc_on','build_olsxc_off','average_rho','build_snxc_on','build_snxc_off','assemble_snxc_on','assemble_snxc_off',
        'build_ca_snxc_on','build_ca_olsxc_on','assemble_h','assemble_mcweda','assemble_hxc','assemble_eh','getenergy','getenergy_hxc','getenergy_mcweda',
        'getenergy_eh','assemble_h_ks','getenergy_KS','assemble_S','assemble_2c_S','assemble_hartree','assemble_scissor','assemble_qmmm',
        'assemble_ca_2c_dip','assemble_ca_3c_dip','assemble_lr_dip','assemble_zw_1c_na','assemble_zw_2c_ct','assemble_zw_3c_ct','assemble_xczw',
        'assemble_zw_off_na','assemble_zw_on_na','build_zw_off_na','build_zw_on_na','assemble_1c_vdip','getenergy_zw']
}

GRID = {
  '' : ['assemble_KS_den0','assemble_KS_den','assemble_KS_usr','laplace_fft','assemble_KS_dcc','assemble_KS_mat','mixer_KS','writeout_charges_KS',
        'writeout_xsf','assemble_KS_vna','get_Hort']
}

INITIALIZERS = {
  '' : ['diagnostics','initatomicE','initconstraints','initcharges','initconstants','initboxes','initkpoints','initmasses','initneighbors',
        'welcome','make_mu2shell','make_munu','make_munuPP','restart','make_munuDipY','make_munuDipX','zero_ang_mom','initamat','make_munuS',
        'initNH','getkpoints','greatKAuto','greatKsubsAuto','initgrid','initdenmat','get_info_orbital','initbasics','initcharges_KS','initcDFT']
}

INTERPOLATERS = {
  '' : ['buildspline_1d','interpolate_1d','interpolate_2d','recover_2c','recover_3c','recover_PP','recoverC','setterp_2d','recover_2cDipY',
        'recover_2cDipX','recover_S','buildspline2_1d','getpsi','getYlm','getvna']
}

LOOPS = {
  '' : ['main_loop','main_loop_MD','main_loop_CG','scf_loop','scf_loop_harris','main_loop_NEB','main_loop_DM','scf_loop_ks','main_loop_importrho',
        'main_loop_TDSE','main_loop_MDET','main_loop_MIN','main_loop_NAC','main_loop_FIRE','main_loop_socket']
}

MAIN = {
  '' : ['fireball']
}

MAIN_SERVER = {
  '' : ['fireball_server']
}

MAIN_SERVER_AMBER = {
  '' : ['fireball_server_amber']
}

MD = {
  '' : ['cross','factorial','predictor','gaussT','corrector','imaged','setgear','phimat','bvec','soldm','NHCThermostat','get_enk','writeHNose',
        'resetNHC','move_ions']
}

MODULES = {
  '' : ['barrier','charges','configuration','constants_fireball','density','dimensions','forces','fragments','gaussG','integrals','interactions',
        'kpoints','neighbor_map','umbrella','steered','optimization','module_dos','dynamo','cproc','noseHoover','scf','grid','wavefunction','neb_module',
        'vnneutral','transport','matmultmod','outputs','options','energy','MD','mpi_main','tdse','bias','nonadiabatic','hartree','sockets','fsockets',
        'fb_socket']
}

MODULES_C = {
  '' : ['classicMD']
}

NEIGHBORS = {
  '' : ['backnay','common_neighbors','find_neigh_max','find_neigh_max_class','mpairnay','neighbors','neighbors_pairs','find_neighPP_max',
        'neighborsPP','common_neighborsPP','num_neigh_tot']
}

PRESSURE = {
  '' : ['hmetric','initpressure','invert3x3']
}

READFILES = {
  '' : ['append_string','read_1c','read_2c','read_3c','readbasis','readbarrier','readdata_2c','readdata_3c','readfragments','readheader_2c',
        'readheader_3c','readinfo','readlvs','readparam','readphi','readpressure','readquench','readsa','readvdw','readcgo','readdos','readgaussG',
        'findFdata','readgrid','read_wf','read_vna','readtrans','readbind','readhop','readdata','readdata_hxc','readdata_mcweda','readdata_xczw',
        'readdata_eh','checksum_options','readdata_KS','getsections','readdata_classicMD','readhartree','readephc']
}

ROTATIONS = {
  '' : ['chooser','chooserd','deps2center','deps3center','makeDmat','makeDmatPP','rotate','rotated','rotatedPP','twister','twisterd',
        'rotatePP','epsilon']
}

THERMOINT = {
  '' : ['cclient']
}

SOCKETS = {
  '' : ['get_geometry','create_socket','send_geometry','sendrecv','soc_init']
}

UMBRELLA = {
  '' : ['assemble_umbrella','Dassemble_umbrella','get_umbrella','readumbrella','assemble_steered','Dassemble_steered','get_steered',
        'readsteered']
}

XC = {
  '' : ['ceperley_alder','cepal']
}

CG = {
  '' : ['cgo','bfgs','l-bfgs-b','FIRE']
}

DOS = {
  '' : ['dos','invierte','writeout_dos','writeout_dosng','hoppings','writeout_atom','hamilt_atom']
}

NEB = {
  '' : ['initneb','neb']
}

TDSE = {
  '' : ['ete_loop','psi2es','eigenHS','tddenmat','tddiag_k','diag_Sk','diag_Hk','allocate_tdse','tdbc','propTpsi','readtdse','initpsi',
        'ortho_H','get_QLow','get_QMul','postete','wrtout_psiT']
}

TRANS = {
  '' : ['assemble_t12_fit','assemble_t12_bare','calcG','assemble_Hsam','assemble_Gsam','assemble_Dxx','sqrt_mat','interpolate_hop',
        'gethop']
}

BIAS = {
  '' : ['assemble_bias','Dassemble_bias','allocate_bias','reallocate_bias','readbias']
}

NAC = {
  '' : ['allocate_nac','assemble_G_S','nacouplings','build_gover1c','init_mdet','mdetdenmat','getforces_mdet','save_mdetstuff','evolve_ks_states',
        'deallocate_nac','dcdt_nac','Dassemble_2c_mdet','Dassemble_2c_PP_mdet','Dassemble_olsxc_on_mdet','Dassemble_olsxc_2c_mdet','Dassemble_3c_mdet',
        'Dassemble_3c_PP_mdet','Dassemble_olsxc_3c_mdet','fewest_switches','mc_switch','transition','Dassemble_ca_2c_mdet','Dassemble_ca_3c_mdet',
        'Dassemble_lr_mdet','Dassemble_ca_olsxc_on_mdet','Dassemble_ca_olsxc_2c_mdet','Dassemble_ca_olsxc_3c_mdet','move_correc','move_predic']
}

QMMM = {
  '' : ['main_loop_MDET_qmmm','main_loop_MD_qmmm','fireball_qmmm_loop']
}

DFTD3 = {
  '' : ['common','sizes','pars','core','api','dftd3_corrections']
}

