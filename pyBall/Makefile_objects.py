
# =================== Files with special compilation 

#INCLUDE = " -I/usr/local/mpich/include "
SPECIAL_CC={
'F77':{ 'lanc', 'blas', 'minilapack', 'qralg', 'invierte' },
'CC' :{ 'sockets', 'cclient' , 'sendrecv', 'soc_init'     },
}

# =================== GROUPS (MODULES) AND OPTIONS (VARIANTS)

all_group_names = ['MODULES','MATH','MAIN','DASSEMBLERS','INTERACTIONS','ALLOCATIONS','ASSEMBLERS','INITIALIZERS','INTERPOLATERS','NEIGHBORS','READFILES','ROTATIONS']

#all_modes   = [ 'DOUBLE', 'GAMMA','OPENMP', 'MPI-k', 'SCALAPACK', 'LAPACK95','ORDERN'   ]
all_modes    = [ 'GAMMA'  ]

optional_modules = []
all_targets = [ 'PROGRAM', 'SERVER']
optionals_modules_dict = dict( optional_modules )

all_optional_modules = [ v for k,vs in optional_modules for v in vs ]
all_variant_names    = [''] + all_modes + all_targets + all_optional_modules



GROUPS = {

'MODULES' : {
'' : [
'dimensions',
'charges',
'configuration',
'constants_fireball',
'density',
'forces',
'integrals',
'interactions',
'kpoints',
'neighbor_map',
'vnneutral',
'wavefunction',
'options',
'energy',
#'scf',
'loops',
'workmat',
'debug',
],
}, #END MODULES

# ======================================
# ========     CORE MODULES
# ======================================

'READFILES' : {
'' : ['append_string','read_1c','read_2c','read_3c','readbasis','readdata_2c','readdata_3c','readheader_2c',
        'readheader_3c','readinfo','readparam',
        'findFdata',
        #'read_wf',
        #'read_vna',
        'readdata','readdata_mcweda',
        'checksum_options',
        'getsections'],
}, #END READFILES

'INTERPOLATERS' : {
'' : ['buildspline_1d','interpolate_1d','interpolate_2d','recover_2c','recover_3c','recover_PP','recoverC','setterp_2d',
#'recover_2cDipY',
#'recover_2cDipX',
'recover_S','buildspline2_1d','getpsi','getYlm','getvna'],
}, #END INTERPOLATERS


'ALLOCATIONS' : {
'' : ['allocate_f','allocate_h','allocate_neigh','allocate_rho','reallocate_f','reallocate_h','reallocate_neigh', 'reallocate_rho'],
}, #END ALLOCATIONS

'INITIALIZERS' : {
'' : [ 
#'diagnostics',
#'initconstraints',
'initcharges',
'initconstants',
'initboxes',
#'initkpoints',
#'initmasses',
'initneighbors',      
#'welcome',
'make_mu2shell',
'make_munu',
'make_munuPP',
#'restart',
#'make_munuDipY',
#'make_munuDipX',
#'zero_ang_mom',
'initamat',
'make_munuS',
#'getkpoints',
#'initdenmat',
'get_info_orbital',
'initbasics'
],
}, #END INITIALIZERS

'NEIGHBORS' : {
'' : ['backnay','common_neighbors','find_neigh_max','mpairnay','neighbors','neighbors_pairs','find_neighPP_max',
        'neighborsPP','common_neighborsPP','num_neigh_tot'],
}, #END NEIGHBORS

#'XC' : {
#'' : ['cepal'],
#}, #END XC


'INTERACTIONS' : {
'':  [
'cl_value', 
'Dtrescentros', 
'doscentros', 
'doscentrosPP', 
#'doscentrosDipY', 
#'doscentrosDipX', 
'get_ewald', 
'trescentros', 
'unocentros', 
'smoother', 
'doscentrosS', 
'trescentrosS', 
'DtrescentrosS', 
'internalLambda', 
'tester2c'
],
}, #END INTERACTIONS

'ASSEMBLERS' : {
'' : [
'assemble_olsxc_1c',
'assemble_2c',
'assemble_3c',
'assemble_ca_2c',
'assemble_3c_PP',
'assemble_2c_PP',
'assemble_ca_3c',
'assemble_F',
'assemble_lr',
'assemble_sVNL',
'assemble_usr',
'buildh',
'assemble_olsxc_on',
'assemble_olsxc_off',
'build_olsxc_on',
'build_olsxc_off',
'average_rho',
'average_ca_rho',
'build_snxc_on',
'build_snxc_off',
'assemble_snxc_on',
'assemble_snxc_off',
'build_ca_snxc_on',
'build_ca_olsxc_on',
'assemble_h',
'assemble_mcweda',
'getenergy',
'getenergy_mcweda',
'assemble_S',
'assemble_2c_S',
#'assemble_ca_2c_dip',
#'assemble_ca_3c_dip',
#'assemble_lr_dip',
#'assemble_1c_vdip',
],
}, #END ASSEMBLERS

'DASSEMBLERS' : {
'' : [
'Dassemble_2c',
'Dassemble_3c',
'Dassemble_ca_2c',
'Dassemble_ca_3c',
'Dassemble_lr',
'Dassemble_snxc_on',
'Dassemble_olsxc_on',
'Dassemble_olsxc_2c',
'Dassemble_olsxc_3c',
'Dassemble_snxc_2c',
'Dassemble_snxc_3c',
'Dassemble_2c_PP',
'Dassemble_3c_PP',
'Dassemble_ca_olsxc_on',
'Dassemble_ca_snxc_on',
'Dassemble_ca_snxc_3c',
'Dassemble_ca_olsxc_3c',
'Dassemble_ca_snxc_2c',
'Dassemble_ca_olsxc_2c',
#'Dassemble_ca_2c_dip',
#'Dassemble_ca_3c_dip',
#'Dassemble_lr_dip',
'getforces_mcweda',
'getforces',
],
}, #END DASSEMBLERS

'MAIN' : {
'' : [
'fireball',
'denmat',
'fermie',
#'build_rho',
#'kspace',
#'kspace2'

'sqrtS',
'ktransform',
'solveH',


],
}, #END FORM_RHO

'MATH' : {
'' : [
'cross',
'factorial',
'cepal'
],
}, #END FORM_RHO

'ROTATIONS' : {
'' : ['chooser','chooserd','deps2center','deps3center','makeDmat','makeDmatPP','rotate','rotated','rotatedPP','twister','twisterd','rotatePP','epsilon'],
}, #END ROTATIONS

# ======================================
# ========     OPTIONAL MODULES
# ======================================

} ######     END GROUPS
