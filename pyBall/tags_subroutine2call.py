





unmatched_calls = []
unmatched_defs  = ['allocate_f', 'allocate_h', 'allocate_neigh', 'assemble_h', 'assemble_hxc_usr', 'assemble_s', 'buildspline2_1d', 'buildspline_1d', 'ceperley_alder', 'dgelements', 'dgelements_vxc', 'dohermit', 'gelements', 'gelements_vxc', 'get_info_orbital', 'getenergy', 'getforces', 'getforces_eh', 'getforces_ks', 'getpsi', 'getvna', 'getylm', 'laplace_fdm', 'setterp_2d', 'test2c']
assemble_1c_vdip = {
'Defs'   : ['ASSEMBLERS/assemble_1c_vdip.f90'],
'Calls'  : ['ASSEMBLERS/assemble_mcweda.f90']
}
assemble_2c = {
'Defs'   : ['ASSEMBLERS/assemble_2c.f90'],
'Calls'  : ['ASSEMBLERS/assemble_eh.f90', 'ASSEMBLERS/assemble_hxc.f90', 'ASSEMBLERS/assemble_mcweda.f90']
}
assemble_2c_pp = {
'Defs'   : ['ASSEMBLERS/assemble_2c_PP.f90'],
'Calls'  : ['ASSEMBLERS/assemble_eh.f90', 'ASSEMBLERS/assemble_hxc.f90', 'ASSEMBLERS/assemble_mcweda.f90']
}
assemble_2c_s = {
'Defs'   : ['ASSEMBLERS/assemble_2c_S.f90'],
'Calls'  : ['ASSEMBLERS/assemble_S.f90']
}
assemble_3c = {
'Defs'   : ['ASSEMBLERS/assemble_3c.f90'],
'Calls'  : ['ASSEMBLERS/assemble_eh.f90', 'ASSEMBLERS/assemble_hxc.f90', 'ASSEMBLERS/assemble_mcweda.f90']
}
assemble_3c_pp = {
'Defs'   : ['ASSEMBLERS/assemble_3c_PP.f90'],
'Calls'  : ['ASSEMBLERS/assemble_eh.f90', 'ASSEMBLERS/assemble_hxc.f90', 'ASSEMBLERS/assemble_mcweda.f90']
}
assemble_ca_2c = {
'Defs'   : ['ASSEMBLERS/assemble_ca_2c.f90'],
'Calls'  : ['ASSEMBLERS/assemble_hxc.f90', 'ASSEMBLERS/assemble_mcweda.f90']
}
assemble_ca_2c_dip = {
'Defs'   : ['ASSEMBLERS/assemble_ca_2c_dip.f90'],
'Calls'  : ['ASSEMBLERS/assemble_hxc.f90', 'ASSEMBLERS/assemble_mcweda.f90']
}
assemble_ca_3c = {
'Defs'   : ['ASSEMBLERS/assemble_ca_3c.f90'],
'Calls'  : ['ASSEMBLERS/assemble_hxc.f90', 'ASSEMBLERS/assemble_mcweda.f90']
}
assemble_ca_3c_dip = {
'Defs'   : ['ASSEMBLERS/assemble_ca_3c_dip.f90'],
'Calls'  : ['ASSEMBLERS/assemble_mcweda.f90']
}
assemble_eh = {
'Defs'   : ['ASSEMBLERS/assemble_eh.f90'],
'Calls'  : ['ASSEMBLERS/assemble_h.f90']
}
assemble_eh_2c = {
'Defs'   : ['ASSEMBLERS/assemble_eh_2c.f90'],
'Calls'  : ['ASSEMBLERS/assemble_eh.f90']
}
assemble_eh_usr = {
'Defs'   : ['ASSEMBLERS/assemble_eh_usr.f90'],
'Calls'  : ['ASSEMBLERS/getenergy_eh.f90']
}
assemble_f = {
'Defs'   : ['ASSEMBLERS/assemble_F.f90'],
'Calls'  : ['DASSEMBLERS/getforces_KS.f90', 'DASSEMBLERS/getforces_eh.f90', 'DASSEMBLERS/getforces_hxc.f90', 'DASSEMBLERS/getforces_mcweda.f90']
}
assemble_hartree = {
'Defs'   : ['ASSEMBLERS/assemble_hartree.f90'],
'Calls'  : ['ASSEMBLERS/assemble_h.f90']
}
assemble_hxc = {
'Defs'   : ['ASSEMBLERS/assemble_hxc.f90'],
'Calls'  : ['ASSEMBLERS/assemble_h.f90']
}
assemble_hxc_1c = {
'Defs'   : ['ASSEMBLERS/assemble_hxc_1c.f90'],
'Calls'  : ['ASSEMBLERS/assemble_eh.f90', 'ASSEMBLERS/assemble_hxc.f90']
}
assemble_hxc_2c = {
'Defs'   : ['ASSEMBLERS/assemble_hxc_2c.f90'],
'Calls'  : ['ASSEMBLERS/assemble_eh.f90', 'ASSEMBLERS/assemble_hxc.f90']
}
assemble_hxc_3c = {
'Defs'   : ['ASSEMBLERS/assemble_hxc_3c.f90'],
'Calls'  : ['ASSEMBLERS/assemble_eh.f90', 'ASSEMBLERS/assemble_hxc.f90']
}
assemble_lr = {
'Defs'   : ['ASSEMBLERS/assemble_lr.f90'],
'Calls'  : ['ASSEMBLERS/assemble_hxc.f90', 'ASSEMBLERS/assemble_mcweda.f90']
}
assemble_lr_dip = {
'Defs'   : ['ASSEMBLERS/assemble_lr_dip.f90'],
'Calls'  : ['ASSEMBLERS/assemble_hxc.f90', 'ASSEMBLERS/assemble_mcweda.f90']
}
assemble_mcweda = {
'Defs'   : ['ASSEMBLERS/assemble_mcweda.f90'],
'Calls'  : ['ASSEMBLERS/assemble_h.f90']
}
assemble_olsxc_1c = {
'Defs'   : ['ASSEMBLERS/assemble_olsxc_1c.f90'],
'Calls'  : ['ASSEMBLERS/assemble_eh.f90', 'ASSEMBLERS/assemble_mcweda.f90']
}
assemble_olsxc_off = {
'Defs'   : ['ASSEMBLERS/assemble_olsxc_off.f90'],
'Calls'  : ['ASSEMBLERS/assemble_eh.f90', 'ASSEMBLERS/assemble_mcweda.f90']
}
assemble_olsxc_on = {
'Defs'   : ['ASSEMBLERS/assemble_olsxc_on.f90'],
'Calls'  : ['ASSEMBLERS/assemble_eh.f90', 'ASSEMBLERS/assemble_mcweda.f90']
}
assemble_scissor = {
'Defs'   : ['ASSEMBLERS/assemble_scissor.f90'],
'Calls'  : ['ASSEMBLERS/assemble_h.f90']
}
assemble_snxc_off = {
'Defs'   : ['ASSEMBLERS/assemble_snxc_off.f90'],
'Calls'  : ['ASSEMBLERS/assemble_eh.f90', 'ASSEMBLERS/assemble_mcweda.f90']
}
assemble_snxc_on = {
'Defs'   : ['ASSEMBLERS/assemble_snxc_on.f90'],
'Calls'  : ['ASSEMBLERS/assemble_eh.f90', 'ASSEMBLERS/assemble_mcweda.f90']
}
assemble_svnl = {
'Defs'   : ['ASSEMBLERS/assemble_sVNL.f90'],
'Calls'  : ['ASSEMBLERS/assemble_eh.f90', 'ASSEMBLERS/assemble_hxc.f90', 'ASSEMBLERS/assemble_mcweda.f90']
}
assemble_usr = {
'Defs'   : ['ASSEMBLERS/assemble_usr.f90'],
'Calls'  : ['ASSEMBLERS/getenergy_eh.f90', 'ASSEMBLERS/getenergy_hxc.f90', 'ASSEMBLERS/getenergy_mcweda.f90']
}
average_ca_rho = {
'Defs'   : ['ASSEMBLERS/average_ca_rho.f90'],
'Calls'  : ['ASSEMBLERS/assemble_mcweda.f90']
}
average_rho = {
'Defs'   : ['ASSEMBLERS/average_rho.f90'],
'Calls'  : ['ASSEMBLERS/assemble_eh.f90', 'ASSEMBLERS/assemble_mcweda.f90']
}
backnay = {
'Defs'   : ['NEIGHBORS/backnay.f90'],
'Calls'  : ['ASSEMBLERS/assemble_S.f90', 'ASSEMBLERS/assemble_eh.f90', 'ASSEMBLERS/assemble_hxc.f90', 'ASSEMBLERS/assemble_mcweda.f90']
}
build_ca_olsxc_on = {
'Defs'   : ['ASSEMBLERS/build_ca_olsxc_on.f90'],
'Calls'  : ['ASSEMBLERS/assemble_olsxc_on.f90']
}
build_ca_snxc_on = {
'Defs'   : ['ASSEMBLERS/build_ca_snxc_on.f90'],
'Calls'  : ['ASSEMBLERS/assemble_snxc_on.f90']
}
build_olsxc_off = {
'Defs'   : ['ASSEMBLERS/build_olsxc_off.f90'],
'Calls'  : ['ASSEMBLERS/assemble_olsxc_off.f90']
}
build_olsxc_on = {
'Defs'   : ['ASSEMBLERS/build_olsxc_on.f90'],
'Calls'  : ['ASSEMBLERS/assemble_olsxc_on.f90']
}
build_snxc_off = {
'Defs'   : ['ASSEMBLERS/build_snxc_off.f90'],
'Calls'  : ['ASSEMBLERS/assemble_snxc_off.f90']
}
build_snxc_on = {
'Defs'   : ['ASSEMBLERS/build_snxc_on.f90'],
'Calls'  : ['ASSEMBLERS/assemble_snxc_on.f90']
}
buildh = {
'Defs'   : ['ASSEMBLERS/buildh.f90'],
'Calls'  : ['ASSEMBLERS/assemble_eh.f90', 'ASSEMBLERS/assemble_hxc.f90', 'ASSEMBLERS/assemble_mcweda.f90']
}
cepal = {
'Defs'   : ['XC/cepal.f90'],
'Calls'  : ['ASSEMBLERS/build_ca_olsxc_on.f90', 'ASSEMBLERS/build_ca_snxc_on.f90', 'ASSEMBLERS/build_olsxc_off.f90', 'ASSEMBLERS/build_olsxc_on.f90', 'ASSEMBLERS/build_snxc_off.f90', 'ASSEMBLERS/build_snxc_on.f90', 'DASSEMBLERS/Dassemble_ca_olsxc_2c.f90', 'DASSEMBLERS/Dassemble_ca_olsxc_3c.f90', 'DASSEMBLERS/Dassemble_ca_olsxc_on.f90', 'DASSEMBLERS/Dassemble_ca_snxc_2c.f90', 'DASSEMBLERS/Dassemble_ca_snxc_3c.f90', 'DASSEMBLERS/Dassemble_ca_snxc_on.f90', 'DASSEMBLERS/Dassemble_olsxc_2c.f90', 'DASSEMBLERS/Dassemble_olsxc_3c.f90', 'DASSEMBLERS/Dassemble_olsxc_on.f90', 'DASSEMBLERS/Dassemble_snxc_2c.f90', 'DASSEMBLERS/Dassemble_snxc_3c.f90', 'DASSEMBLERS/Dassemble_snxc_on.f90']
}
chooser = {
'Defs'   : ['ROTATIONS/chooser.f90'],
'Calls'  : ['ROTATIONS/makeDmat.f90', 'ROTATIONS/makeDmatPP.f90', 'ROTATIONS/rotate.f90', 'ROTATIONS/rotatePP.f90']
}
chooserd = {
'Defs'   : ['ROTATIONS/chooserd.f90'],
'Calls'  : ['ROTATIONS/makeDmat.f90', 'ROTATIONS/makeDmatPP.f90']
}
cl_value = {
'Defs'   : ['INTERACTIONS/cl_value.f90'],
'Calls'  : ['ASSEMBLERS/assemble_2c_PP.f90', 'ASSEMBLERS/assemble_3c_PP.f90', 'DASSEMBLERS/Dassemble_2c_PP.f90', 'DASSEMBLERS/Dassemble_3c_PP.f90']
}
common_neighbors = {
'Defs'   : ['NEIGHBORS/common_neighbors.f90'],
'Calls'  : ['ASSEMBLERS/assemble_S.f90', 'ASSEMBLERS/assemble_eh.f90', 'ASSEMBLERS/assemble_hxc.f90', 'ASSEMBLERS/assemble_mcweda.f90']
}
common_neighborspp = {
'Defs'   : ['NEIGHBORS/common_neighborsPP.f90'],
'Calls'  : ['ASSEMBLERS/assemble_eh.f90', 'ASSEMBLERS/assemble_hxc.f90', 'ASSEMBLERS/assemble_mcweda.f90']
}
cross = {
'Defs'   : ['MATH/cross.f90'],
'Calls'  : ['INTERACTIONS/get_ewald.f90', 'INTERACTIONS/get_ewald_OMP.f90']
}
dassemble_2c = {
'Defs'   : ['DASSEMBLERS/Dassemble_2c.f90'],
'Calls'  : ['DASSEMBLERS/getforces_KS.f90', 'DASSEMBLERS/getforces_eh.f90', 'DASSEMBLERS/getforces_hxc.f90', 'DASSEMBLERS/getforces_mcweda.f90']
}
dassemble_2c_pp = {
'Defs'   : ['DASSEMBLERS/Dassemble_2c_PP.f90'],
'Calls'  : ['DASSEMBLERS/getforces_KS.f90', 'DASSEMBLERS/getforces_eh.f90', 'DASSEMBLERS/getforces_hxc.f90', 'DASSEMBLERS/getforces_mcweda.f90']
}
dassemble_3c = {
'Defs'   : ['DASSEMBLERS/Dassemble_3c.f90'],
'Calls'  : ['DASSEMBLERS/getforces_KS.f90', 'DASSEMBLERS/getforces_eh.f90', 'DASSEMBLERS/getforces_hxc.f90', 'DASSEMBLERS/getforces_mcweda.f90']
}
dassemble_3c_pp = {
'Defs'   : ['DASSEMBLERS/Dassemble_3c_PP.f90'],
'Calls'  : ['DASSEMBLERS/getforces_KS.f90', 'DASSEMBLERS/getforces_eh.f90', 'DASSEMBLERS/getforces_hxc.f90', 'DASSEMBLERS/getforces_mcweda.f90']
}
dassemble_ca_2c = {
'Defs'   : ['DASSEMBLERS/Dassemble_ca_2c.f90'],
'Calls'  : ['DASSEMBLERS/getforces_hxc.f90', 'DASSEMBLERS/getforces_mcweda.f90']
}
dassemble_ca_2c_dip = {
'Defs'   : ['DASSEMBLERS/Dassemble_ca_2c_dip.f90'],
'Calls'  : ['DASSEMBLERS/getforces_hxc.f90', 'DASSEMBLERS/getforces_mcweda.f90']
}
dassemble_ca_3c = {
'Defs'   : ['DASSEMBLERS/Dassemble_ca_3c.f90'],
'Calls'  : ['DASSEMBLERS/getforces_hxc.f90', 'DASSEMBLERS/getforces_mcweda.f90']
}
dassemble_ca_3c_dip = {
'Defs'   : ['DASSEMBLERS/Dassemble_ca_3c_dip.f90'],
'Calls'  : ['DASSEMBLERS/getforces_hxc.f90', 'DASSEMBLERS/getforces_mcweda.f90']
}
dassemble_ca_olsxc_2c = {
'Defs'   : ['DASSEMBLERS/Dassemble_ca_olsxc_2c.f90'],
'Calls'  : ['DASSEMBLERS/getforces_mcweda.f90']
}
dassemble_ca_olsxc_3c = {
'Defs'   : ['DASSEMBLERS/Dassemble_ca_olsxc_3c.f90'],
'Calls'  : ['DASSEMBLERS/getforces_mcweda.f90']
}
dassemble_ca_olsxc_on = {
'Defs'   : ['DASSEMBLERS/Dassemble_ca_olsxc_on.f90'],
'Calls'  : ['DASSEMBLERS/getforces_mcweda.f90']
}
dassemble_ca_snxc_2c = {
'Defs'   : ['DASSEMBLERS/Dassemble_ca_snxc_2c.f90'],
'Calls'  : ['DASSEMBLERS/getforces_mcweda.f90']
}
dassemble_ca_snxc_3c = {
'Defs'   : ['DASSEMBLERS/Dassemble_ca_snxc_3c.f90'],
'Calls'  : ['DASSEMBLERS/getforces_mcweda.f90']
}
dassemble_ca_snxc_on = {
'Defs'   : ['DASSEMBLERS/Dassemble_ca_snxc_on.f90'],
'Calls'  : ['DASSEMBLERS/getforces_mcweda.f90']
}
dassemble_eh_2c = {
'Defs'   : ['DASSEMBLERS/Dassemble_eh_2c.f90'],
'Calls'  : ['DASSEMBLERS/getforces_eh.f90']
}
dassemble_hxc_2c = {
'Defs'   : ['DASSEMBLERS/Dassemble_hxc_2c.f90'],
'Calls'  : ['DASSEMBLERS/getforces_hxc.f90']
}
dassemble_hxc_3c = {
'Defs'   : ['DASSEMBLERS/Dassemble_hxc_3c.f90'],
'Calls'  : ['DASSEMBLERS/getforces_hxc.f90']
}
dassemble_lr = {
'Defs'   : ['DASSEMBLERS/Dassemble_lr.f90', 'DASSEMBLERS/Dassemble_lr_OMP.f90'],
'Calls'  : ['DASSEMBLERS/getforces_hxc.f90', 'DASSEMBLERS/getforces_mcweda.f90']
}
dassemble_lr_dip = {
'Defs'   : ['DASSEMBLERS/Dassemble_lr_dip.f90'],
'Calls'  : ['DASSEMBLERS/getforces_hxc.f90', 'DASSEMBLERS/getforces_mcweda.f90']
}
dassemble_olsxc_2c = {
'Defs'   : ['DASSEMBLERS/Dassemble_olsxc_2c.f90'],
'Calls'  : ['DASSEMBLERS/getforces_eh.f90', 'DASSEMBLERS/getforces_mcweda.f90']
}
dassemble_olsxc_3c = {
'Defs'   : ['DASSEMBLERS/Dassemble_olsxc_3c.f90'],
'Calls'  : ['DASSEMBLERS/getforces_eh.f90', 'DASSEMBLERS/getforces_mcweda.f90']
}
dassemble_olsxc_on = {
'Defs'   : ['DASSEMBLERS/Dassemble_olsxc_on.f90'],
'Calls'  : ['DASSEMBLERS/getforces_eh.f90', 'DASSEMBLERS/getforces_mcweda.f90']
}
dassemble_snxc_2c = {
'Defs'   : ['DASSEMBLERS/Dassemble_snxc_2c.f90'],
'Calls'  : ['DASSEMBLERS/getforces_eh.f90', 'DASSEMBLERS/getforces_mcweda.f90']
}
dassemble_snxc_3c = {
'Defs'   : ['DASSEMBLERS/Dassemble_snxc_3c.f90'],
'Calls'  : ['DASSEMBLERS/getforces_eh.f90', 'DASSEMBLERS/getforces_mcweda.f90']
}
dassemble_snxc_on = {
'Defs'   : ['DASSEMBLERS/Dassemble_snxc_on.f90'],
'Calls'  : ['DASSEMBLERS/getforces_eh.f90', 'DASSEMBLERS/getforces_mcweda.f90']
}
deps2cent = {
'Defs'   : ['ROTATIONS/deps2center.f90'],
'Calls'  : ['ASSEMBLERS/assemble_2c.f90', 'ASSEMBLERS/assemble_2c_S.f90', 'ASSEMBLERS/assemble_ca_2c.f90', 'ASSEMBLERS/assemble_ca_2c_dip.f90', 'ASSEMBLERS/assemble_ca_3c.f90', 'ASSEMBLERS/assemble_ca_3c_dip.f90', 'ASSEMBLERS/assemble_hxc_2c.f90', 'ASSEMBLERS/assemble_olsxc_off.f90', 'ASSEMBLERS/assemble_sVNL.f90', 'ASSEMBLERS/average_ca_rho.f90', 'ASSEMBLERS/average_rho.f90', 'DASSEMBLERS/Dassemble_2c.f90', 'DASSEMBLERS/Dassemble_ca_2c.f90', 'DASSEMBLERS/Dassemble_ca_2c_dip.f90', 'DASSEMBLERS/Dassemble_ca_olsxc_2c.f90', 'DASSEMBLERS/Dassemble_ca_snxc_2c.f90', 'DASSEMBLERS/Dassemble_hxc_2c.f90', 'DASSEMBLERS/Dassemble_olsxc_2c.f90', 'DASSEMBLERS/Dassemble_snxc_2c.f90']
}
deps3center = {
'Defs'   : ['ROTATIONS/deps3center.f90'],
'Calls'  : ['DASSEMBLERS/Dassemble_3c.f90', 'DASSEMBLERS/Dassemble_ca_3c.f90', 'DASSEMBLERS/Dassemble_ca_3c_dip.f90', 'DASSEMBLERS/Dassemble_ca_olsxc_3c.f90', 'DASSEMBLERS/Dassemble_ca_snxc_3c.f90', 'DASSEMBLERS/Dassemble_hxc_3c.f90', 'DASSEMBLERS/Dassemble_olsxc_3c.f90', 'DASSEMBLERS/Dassemble_snxc_3c.f90']
}
doscentros = {
'Defs'   : ['INTERACTIONS/doscentros.f90'],
'Calls'  : ['ASSEMBLERS/assemble_2c.f90', 'ASSEMBLERS/assemble_2c_S.f90', 'ASSEMBLERS/assemble_ca_2c.f90', 'ASSEMBLERS/assemble_ca_2c_dip.f90', 'ASSEMBLERS/assemble_hxc_2c.f90', 'ASSEMBLERS/assemble_olsxc_off.f90', 'ASSEMBLERS/average_ca_rho.f90', 'ASSEMBLERS/average_rho.f90', 'DASSEMBLERS/Dassemble_2c.f90', 'DASSEMBLERS/Dassemble_ca_2c.f90', 'DASSEMBLERS/Dassemble_ca_2c_dip.f90', 'DASSEMBLERS/Dassemble_ca_olsxc_2c.f90', 'DASSEMBLERS/Dassemble_hxc_2c.f90', 'DASSEMBLERS/Dassemble_olsxc_2c.f90']
}
doscentrosdipx = {
'Defs'   : ['INTERACTIONS/doscentrosDipX.f90'],
'Calls'  : ['ASSEMBLERS/assemble_2c.f90']
}
doscentrosdipy = {
'Defs'   : ['INTERACTIONS/doscentrosDipY.f90'],
'Calls'  : ['ASSEMBLERS/assemble_2c.f90']
}
doscentrospp = {
'Defs'   : ['INTERACTIONS/doscentrosPP.f90'],
'Calls'  : ['ASSEMBLERS/assemble_sVNL.f90']
}
doscentross = {
'Defs'   : ['INTERACTIONS/doscentrosS.f90'],
'Calls'  : ['ASSEMBLERS/average_ca_rho.f90', 'ASSEMBLERS/average_rho.f90']
}
dtrescentros = {
'Defs'   : ['INTERACTIONS/Dtrescentros.f90'],
'Calls'  : ['DASSEMBLERS/Dassemble_3c.f90', 'DASSEMBLERS/Dassemble_ca_3c.f90', 'DASSEMBLERS/Dassemble_ca_3c_dip.f90', 'DASSEMBLERS/Dassemble_ca_olsxc_3c.f90', 'DASSEMBLERS/Dassemble_ca_snxc_3c.f90', 'DASSEMBLERS/Dassemble_hxc_3c.f90', 'DASSEMBLERS/Dassemble_olsxc_3c.f90', 'DASSEMBLERS/Dassemble_snxc_3c.f90']
}
dtrescentross = {
'Defs'   : ['INTERACTIONS/DtrescentrosS.f90'],
'Calls'  : ['DASSEMBLERS/Dassemble_ca_olsxc_3c.f90', 'DASSEMBLERS/Dassemble_ca_snxc_3c.f90', 'DASSEMBLERS/Dassemble_olsxc_3c.f90', 'DASSEMBLERS/Dassemble_snxc_3c.f90']
}
epsilon = {
'Defs'   : ['ROTATIONS/epsilon.f90'],
'Calls'  : ['ASSEMBLERS/assemble_2c.f90', 'ASSEMBLERS/assemble_2c_S.f90', 'ASSEMBLERS/assemble_3c.f90', 'ASSEMBLERS/assemble_ca_2c.f90', 'ASSEMBLERS/assemble_ca_2c_dip.f90', 'ASSEMBLERS/assemble_ca_3c.f90', 'ASSEMBLERS/assemble_ca_3c_dip.f90', 'ASSEMBLERS/assemble_hxc_2c.f90', 'ASSEMBLERS/assemble_hxc_3c.f90', 'ASSEMBLERS/assemble_olsxc_off.f90', 'ASSEMBLERS/assemble_sVNL.f90', 'ASSEMBLERS/average_ca_rho.f90', 'ASSEMBLERS/average_rho.f90', 'DASSEMBLERS/Dassemble_2c.f90', 'DASSEMBLERS/Dassemble_3c.f90', 'DASSEMBLERS/Dassemble_ca_2c.f90', 'DASSEMBLERS/Dassemble_ca_2c_dip.f90', 'DASSEMBLERS/Dassemble_ca_3c.f90', 'DASSEMBLERS/Dassemble_ca_3c_dip.f90', 'DASSEMBLERS/Dassemble_ca_olsxc_2c.f90', 'DASSEMBLERS/Dassemble_ca_olsxc_3c.f90', 'DASSEMBLERS/Dassemble_ca_snxc_2c.f90', 'DASSEMBLERS/Dassemble_ca_snxc_3c.f90', 'DASSEMBLERS/Dassemble_hxc_2c.f90', 'DASSEMBLERS/Dassemble_hxc_3c.f90', 'DASSEMBLERS/Dassemble_olsxc_2c.f90', 'DASSEMBLERS/Dassemble_olsxc_3c.f90', 'DASSEMBLERS/Dassemble_snxc_2c.f90', 'DASSEMBLERS/Dassemble_snxc_3c.f90']
}
fillneigh = {
'Defs'   : ['NEIGHBORS/neighbors.f90'],
'Calls'  : ['NEIGHBORS/neighbors.f90']
}
fillneigh_class = {
'Defs'   : ['NEIGHBORS/neighbors.f90'],
'Calls'  : ['NEIGHBORS/neighbors.f90']
}
find_neigh_max = {
'Defs'   : ['NEIGHBORS/find_neigh_max.f90'],
'Calls'  : ['ALLOCATIONS/allocate_neigh.f90', 'ALLOCATIONS/reallocate_neigh.f90']
}
find_neighpp_max = {
'Defs'   : ['NEIGHBORS/find_neighPP_max.f90'],
'Calls'  : ['ALLOCATIONS/allocate_neigh.f90', 'ALLOCATIONS/reallocate_neigh.f90']
}
get_atom_indices = {
'Defs'   : ['INTERACTIONS/get_ewald.f90', 'INTERACTIONS/get_ewald_OMP.f90'],
'Calls'  : ['INTERACTIONS/get_ewald.f90', 'INTERACTIONS/get_ewald_OMP.f90']
}
get_ewald = {
'Defs'   : ['INTERACTIONS/get_ewald.f90', 'INTERACTIONS/get_ewald_OMP.f90'],
'Calls'  : ['ASSEMBLERS/assemble_eh.f90', 'ASSEMBLERS/assemble_hxc.f90', 'ASSEMBLERS/assemble_mcweda.f90', 'ASSEMBLERS/getenergy_eh.f90', 'ASSEMBLERS/getenergy_hxc.f90', 'ASSEMBLERS/getenergy_mcweda.f90']
}
get_vdw = {
'Defs'   : ['INTERACTIONS/get_vdw.f90'],
'Calls'  : ['ASSEMBLERS/getenergy_eh.f90', 'ASSEMBLERS/getenergy_hxc.f90', 'ASSEMBLERS/getenergy_mcweda.f90']
}
getenergy_eh = {
'Defs'   : ['ASSEMBLERS/getenergy_eh.f90'],
'Calls'  : ['ASSEMBLERS/getenergy.f90']
}
getenergy_hxc = {
'Defs'   : ['ASSEMBLERS/getenergy_hxc.f90'],
'Calls'  : ['ASSEMBLERS/getenergy.f90']
}
getenergy_mcweda = {
'Defs'   : ['ASSEMBLERS/getenergy_mcweda.f90'],
'Calls'  : ['ASSEMBLERS/getenergy.f90']
}
getforces_hxc = {
'Defs'   : ['DASSEMBLERS/getforces_hxc.f90'],
'Calls'  : ['DASSEMBLERS/getforces.f90']
}
getforces_mcweda = {
'Defs'   : ['DASSEMBLERS/getforces_mcweda.f90'],
'Calls'  : ['DASSEMBLERS/getforces.f90']
}
getharmonic = {
'Defs'   : ['INTERACTIONS/getHarmonic.f90'],
'Calls'  : ['ASSEMBLERS/getenergy_eh.f90', 'ASSEMBLERS/getenergy_hxc.f90', 'ASSEMBLERS/getenergy_mcweda.f90']
}
initneighbors = {
'Defs'   : ['INITIALIZERS/initneighbors.f90'],
'Calls'  : ['ASSEMBLERS/assemble_S.f90', 'ASSEMBLERS/assemble_eh.f90', 'ASSEMBLERS/assemble_hxc.f90', 'ASSEMBLERS/assemble_mcweda.f90']
}
interpolate_1d = {
'Defs'   : ['INTERPOLATERS/interpolate_1d.f90'],
'Calls'  : ['ASSEMBLERS/assemble_eh_2c.f90', 'ASSEMBLERS/assemble_eh_usr.f90', 'ASSEMBLERS/assemble_hxc_usr.f90', 'ASSEMBLERS/assemble_usr.f90', 'DASSEMBLERS/Dassemble_eh_2c.f90', 'INTERACTIONS/doscentros.f90', 'INTERACTIONS/doscentrosDipX.f90', 'INTERACTIONS/doscentrosDipY.f90', 'INTERACTIONS/doscentrosPP.f90', 'INTERACTIONS/doscentrosS.f90', 'INTERACTIONS/tester2c.f90']
}
interpolate_2d = {
'Defs'   : ['INTERPOLATERS/interpolate_2d.f90'],
'Calls'  : ['INTERACTIONS/Dtrescentros.f90', 'INTERACTIONS/DtrescentrosS.f90', 'INTERACTIONS/trescentros.f90', 'INTERACTIONS/trescentrosS.f90']
}
makedmat = {
'Defs'   : ['ROTATIONS/makeDmat.f90'],
'Calls'  : ['ROTATIONS/rotated.f90']
}
makedmatpp = {
'Defs'   : ['ROTATIONS/makeDmatPP.f90'],
'Calls'  : ['ROTATIONS/rotatedPP.f90']
}
neighbors = {
'Defs'   : ['NEIGHBORS/neighbors.f90'],
'Calls'  : ['ASSEMBLERS/assemble_S.f90', 'ASSEMBLERS/assemble_eh.f90', 'ASSEMBLERS/assemble_hxc.f90', 'ASSEMBLERS/assemble_mcweda.f90']
}
neighbors_pairs = {
'Defs'   : ['NEIGHBORS/neighbors_pairs.f90'],
'Calls'  : ['ASSEMBLERS/assemble_S.f90', 'ASSEMBLERS/assemble_eh.f90', 'ASSEMBLERS/assemble_hxc.f90', 'ASSEMBLERS/assemble_mcweda.f90']
}
neighborspp = {
'Defs'   : ['NEIGHBORS/neighborsPP.f90'],
'Calls'  : ['ASSEMBLERS/assemble_eh.f90', 'ASSEMBLERS/assemble_hxc.f90', 'ASSEMBLERS/assemble_mcweda.f90']
}
num_neigh_tot = {
'Defs'   : ['NEIGHBORS/num_neigh_tot.f90'],
'Calls'  : ['ASSEMBLERS/assemble_S.f90', 'ASSEMBLERS/assemble_eh.f90', 'ASSEMBLERS/assemble_hxc.f90', 'ASSEMBLERS/assemble_mcweda.f90']
}
reallocate_2d_iarray = {
'Defs'   : ['NEIGHBORS/neighbors.f90'],
'Calls'  : ['NEIGHBORS/neighbors.f90']
}
reallocate_f = {
'Defs'   : ['ALLOCATIONS/reallocate_f.f90'],
'Calls'  : ['ALLOCATIONS/reallocate_neigh.f90']
}
reallocate_h = {
'Defs'   : ['ALLOCATIONS/reallocate_h.f90'],
'Calls'  : ['ALLOCATIONS/reallocate_neigh.f90']
}
reallocate_neigh = {
'Defs'   : ['ALLOCATIONS/reallocate_neigh.f90'],
'Calls'  : ['ASSEMBLERS/assemble_S.f90', 'ASSEMBLERS/assemble_eh.f90', 'ASSEMBLERS/assemble_hxc.f90', 'ASSEMBLERS/assemble_mcweda.f90']
}
reallocate_rho = {
'Defs'   : ['ALLOCATIONS/reallocate_rho.f90'],
'Calls'  : ['ALLOCATIONS/reallocate_neigh.f90']
}
recover_2c = {
'Defs'   : ['INTERPOLATERS/recover_2c.f90'],
'Calls'  : ['INTERACTIONS/doscentros.f90']
}
recover_2cdipx = {
'Defs'   : ['INTERPOLATERS/recover_2cDipX.f90'],
'Calls'  : ['INTERACTIONS/doscentrosDipX.f90']
}
recover_2cdipy = {
'Defs'   : ['INTERPOLATERS/recover_2cDipY.f90'],
'Calls'  : ['INTERACTIONS/doscentrosDipY.f90']
}
recover_3c = {
'Defs'   : ['INTERPOLATERS/recover_3c.f90'],
'Calls'  : ['INTERACTIONS/Dtrescentros.f90', 'INTERACTIONS/trescentros.f90']
}
recover_pp = {
'Defs'   : ['INTERPOLATERS/recover_PP.f90'],
'Calls'  : ['INTERACTIONS/doscentrosPP.f90']
}
recover_s = {
'Defs'   : ['INTERPOLATERS/recover_S.f90'],
'Calls'  : ['INTERACTIONS/DtrescentrosS.f90', 'INTERACTIONS/doscentrosS.f90', 'INTERACTIONS/trescentrosS.f90']
}
recoverc = {
'Defs'   : ['INTERPOLATERS/recoverC.f90'],
'Calls'  : ['ASSEMBLERS/assemble_eh_2c.f90', 'ASSEMBLERS/assemble_eh_usr.f90', 'ASSEMBLERS/assemble_usr.f90', 'DASSEMBLERS/Dassemble_eh_2c.f90']
}
rotate_fb = {
'Defs'   : ['ROTATIONS/rotate.f90'],
'Calls'  : ['INTERACTIONS/Dtrescentros.f90', 'INTERACTIONS/doscentros.f90', 'INTERACTIONS/doscentrosDipX.f90', 'INTERACTIONS/doscentrosDipY.f90', 'INTERACTIONS/trescentros.f90']
}
rotated = {
'Defs'   : ['ROTATIONS/rotated.f90'],
'Calls'  : ['INTERACTIONS/Dtrescentros.f90', 'INTERACTIONS/doscentros.f90', 'INTERACTIONS/doscentrosDipX.f90', 'INTERACTIONS/doscentrosDipY.f90']
}
rotatedpp = {
'Defs'   : ['ROTATIONS/rotatedPP.f90'],
'Calls'  : ['INTERACTIONS/doscentrosPP.f90']
}
rotatepp = {
'Defs'   : ['ROTATIONS/rotatePP.f90'],
'Calls'  : ['INTERACTIONS/doscentrosPP.f90']
}
smoother = {
'Defs'   : ['INTERACTIONS/smoother.f90'],
'Calls'  : ['ASSEMBLERS/assemble_ca_2c.f90', 'ASSEMBLERS/assemble_ca_2c_dip.f90', 'ASSEMBLERS/assemble_ca_3c.f90', 'ASSEMBLERS/assemble_ca_3c_dip.f90', 'DASSEMBLERS/Dassemble_ca_2c.f90', 'DASSEMBLERS/Dassemble_ca_2c_dip.f90', 'DASSEMBLERS/Dassemble_ca_3c.f90', 'DASSEMBLERS/Dassemble_ca_3c_dip.f90']
}
trescentros = {
'Defs'   : ['INTERACTIONS/trescentros.f90'],
'Calls'  : ['ASSEMBLERS/assemble_hxc_3c.f90', 'ASSEMBLERS/average_ca_rho.f90', 'ASSEMBLERS/average_rho.f90']
}
trescentross = {
'Defs'   : ['INTERACTIONS/trescentrosS.f90'],
'Calls'  : ['ASSEMBLERS/average_ca_rho.f90', 'ASSEMBLERS/average_rho.f90', 'DASSEMBLERS/Dassemble_ca_snxc_3c.f90', 'DASSEMBLERS/Dassemble_olsxc_3c.f90', 'DASSEMBLERS/Dassemble_snxc_3c.f90']
}
twister = {
'Defs'   : ['ROTATIONS/twister.f90'],
'Calls'  : ['ROTATIONS/rotate.f90', 'ROTATIONS/rotatePP.f90', 'ROTATIONS/rotated.f90', 'ROTATIONS/rotatedPP.f90']
}
twisterd = {
'Defs'   : ['ROTATIONS/twisterd.f90'],
'Calls'  : ['ROTATIONS/rotated.f90', 'ROTATIONS/rotatedPP.f90']
}
unocentros = {
'Defs'   : ['INTERACTIONS/unocentros.f90'],
'Calls'  : ['ASSEMBLERS/assemble_hxc_1c.f90', 'ASSEMBLERS/assemble_olsxc_1c.f90']
}
