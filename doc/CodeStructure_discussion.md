
## Lets define minimal necessary core:

I think the ***services*** around like `READFILES`, `LOOPS`,  `SOLVESH`, `FORM_RHO`,  `INITIALIZERS`, `ALLOCATIONS`, `NEIGHBORS`, `MODULES`  etc. are secondary … most of this can be better replaced by some functions from numpy/scipy/ASE  … the rest can be easily adjusted after we define the ***core***. 

What I mean by the ***core*** is the minimal set of functions you need to **assemble Fireball Hamiltonian** (because this is the part special for Fireball, which cannot be found in some general numerical package). (And then corresponding forces, since I want relax molecular geometry).  And there I’m a bit lost, because I don’t know what each of these assemblers does, and which are needed always and which are supplementary ( also it depends on settings like `idogs` and `iharris` `mcWEDA` `iTheory` etc.  …. we should choose just one for the start (can add other options later) )

So if somebody can say which of the files are ***essential***  and which are ***optional/supplementary***?


- `ASSEMBLERS`  / `assemble_*.f90`  …  These are for Hamiltonian
- `DASSEMBLERS` /`Dassemble_*.f90` … These are for derivative of Hamiltonian => Forces … (Right?). Since each `Dassemble_*.f90` is connected to particular `assemble_*.f90` we don’t have to discuss them separately
- `INTERACTIONS` … These are functions called from `ASSEMBLERS` and `DASSEMBLERS` for each matrix element? Right? Are they all essential ?
- `ROTATIONS`    … This is probably all essential (?)  … actually this is one thing which I would like to reuse for many other things (perhaps rewrite it to C/OpenCL)
- `INTERPOLATERS` … perhaps need them all (?) to read form `Fdata` tables.
# There are shortcuts, what they mean:

(also if you can point out if they are general, or special for some theory-setting)
_1c, _2c, _3c  … 1 centre, 2 center, 3 center
_ca_2c, _ca_3c  …. ?
_vdip, _dip   … something about dipole (?)
_eh …. Hartree energy?
_usr    …. ?
_F      … ?
_Ir     … ?
mcweda … mcWeda
KS … KhonSham
_S … overlap matrix
_olsxc,  _snxc, _xczw    …. different kinds of exchange correlation functionals ?
_zw  …. ?
_off, _on  … offsite, onsite (?)
_VNA … neutral atom potential
_VXC  … exchange correlation potential
_VNL  … non-local potential?

# PLEASE SELECT WHICH CAN BE REMOVED:


## ASSEMBLERS


[ ] assemble_1c_vdip.f90
[ ] assemble_2c.f90
[x] assemble_2c_ordern_final.f90
[x] assemble_2c_ordern_init.f90
[ ] assemble_2c_PP.f90
[ ] assemble_2c_S.f90
[ ] assemble_3c.f90
[x] assemble_3c_ordern_final.f90
[ ] assemble_3c_PP.f90
[ ] assemble_ca_2c_dip.f90
[ ] assemble_ca_2c.f90
[x] assemble_ca_2c_ordern_final.f90
[ ] assemble_ca_3c_dip.f90
[ ] assemble_ca_3c.f90
[x] assemble_ca_3c_ordern_final.f90
[ ] assemble_eh_2c.f90
[x] assemble_eh_2c_ordern_final.f90
[ ] assemble_eh.f90
[ ] assemble_eh_usr.f90
[ ] assemble_F.f90
[ ] assemble_hartree.f90
[ ] assemble_h.f90
[ ] assemble_h_ks.f90
[ ] assemble_hxc_1c.f90
[ ] assemble_hxc_2c.f90
[ ] assemble_hxc_3c.f90
[ ] assemble_hxc.f90
[ ] assemble_hxc_usr.f90
[ ] assemble_lr_dip.f90
[ ] assemble_lr.f90
[x] assemble_lr_ordern_final.f90
[ ] assemble_mcweda.f90
[ ] assemble_olsxc_1c.f90
[ ] assemble_olsxc_off.f90
[ ] assemble_olsxc_on.f90
[x] assemble_scissor.f90
[ ] assemble_S.f90
[ ] assemble_snxc_off.f90
[ ] assemble_snxc_on.f90
[ ] assemble_sVNL.f90
[ ] assemble_usr.f90
[ ] assemble_xczw.f90
[ ] assemble_zw_1c_na.f90
[ ] assemble_zw_2c_ct.f90
[ ] assemble_zw_3c_ct.f90
[ ] assemble_zw_off_na.f90
[ ] assemble_zw_on_na.f90
[ ] average_ca_rho.f90
[ ] average_rho.f90
[ ] build_ca_olsxc_on.f90
[ ] build_ca_snxc_on.f90
[ ] buildh.f90
[ ] build_olsxc_off.f90
[ ] build_olsxc_on.f90
[ ] build_snxc_off.f90
[ ] build_snxc_on.f90
[ ] build_zw_off_na.f90
[ ] build_zw_on_na.f90
[ ] getenergy_eh.f90
[ ] getenergy.f90
[ ] getenergy_hxc.f90
[ ] getenergy_KS.f90
[ ] getenergy_mcweda.f90
[ ] getenergy_zw.f90
[ ] getherm.f90
[ ] laplace_fdm.f90


## INTERACTIONS


[ ] cl_value.f90
[ ] Dgelements.f90
[ ] DgelementsG_overlap.f90
[ ] DgelementsGS_overlap.f90
[ ] DgelementsGS_VXC.f90
[ ] DgelementsG_VNA.f90
[ ] DgelementsG_VNA_SH.f90
[ ] DgelementsG_VNA_SH_RNA.f90
[ ] DgelementsG_VXC.f90
[ ] Dgelements_VXC.f90
[ ] doscentrosDipX.f90
[ ] doscentrosDipY.f90
[ ] doscentros.f90
[ ] doscentrosG_overlap.f90
[ ] doscentrosGS_overlap.f90
[ ] doscentrosPP.f90
[ ] doscentrosS.f90
[ ] dosgaussians.f90
[ ] Dtrescentros.f90
[ ] DtrescentrosGHXC_VXC.f90
[ ] DtrescentrosGS_VXC.f90
[ ] DtrescentrosG_VNA.f90
[ ] DtrescentrosG_VNA_SH.f90
[ ] DtrescentrosG_VXC.f90
[ ] DtrescentrosS.f90
[ ] gelements.f90
[ ] gelementsG_overlap.f90
[ ] gelementsGS_overlap.f90
[ ] gelementsGS_VXC.f90
[ ] gelementsG_VNA.f90
[ ] gelementsG_VNA_SH.f90
[ ] gelementsG_VNA_SH_RNA.f90
[ ] gelementsG_VXC.f90
[ ] gelements_VXC.f90
[ ] get_ewald.f90
[ ] get_ewald_OMP.f90
[ ] getHarmonic.f90
[ ] get_info_orbital.f90
[ ] get_vdw.f90
[ ] internalLambda.f90
[ ] smoother.f90
[ ] tester2c.f90
[ ] trescentros.f90
[ ] trescentrosGHXC_VXC.f90
[ ] trescentrosGS_VXC.f90
[ ] trescentrosG_VNA.f90
[ ] trescentrosG_VNA_SH.f90
[ ] trescentrosG_VXC.f90
[ ] trescentrosS.f90
[ ] unocentros.f90


## ROTATIONS


[ ] chooserd.f90
[ ] chooser.f90
[ ] deps2center.f90
[ ] deps3center.f90
[ ] epsilon.f90
[ ] makeDmat.f90
[ ] makeDmatPP.f90
[ ] rotated.f90
[ ] rotatedPP.f90
[ ] rotate.f90
[ ] rotatePP.f90
[ ] twisterd.f90
[ ] twister.f90


## INTERPOLATERS


[ ] buildspline_1d.f90
[ ] buildspline2_1d.f90
[ ] getpsi.f90
[ ] getvna.f90
[ ] getYlm.f90
[ ] interpolate_1d.f90
[ ] interpolate_2d.f90
[ ] recover_2cDipX.f90
[ ] recover_2cDipY.f90
[ ] recover_2c.f90
[ ] recover_3c.f90
[ ] recoverC.f90
[ ] recover_PP.f90
[ ] recover_S.f90
[ ] setterp_2d.f90



