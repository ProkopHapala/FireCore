INCLUDES = 
LFLAGS   = -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -Bdynamic

F90 = gfortran
FFLAGS =   -fPIC   -Ofast -march=native -mtune=native -fPIC -freal-4-real-8  -ffree-form -ffree-line-length-none 
FFLAGS77 =   -fPIC   -Ofast -march=native -mtune=native -fPIC -freal-4-real-8  
FFLAGSC =   -fPIC   -Ofast -march=native -mtune=native -fPIC -freal-4-real-8  -ffree-form -ffree-line-length-none 
MODULES =  dimensions.o charges.o configuration.o constants_fireball.o density.o forces.o Fdata3c.o integrals.o interactions.o \
	 kpoints.o neighbor_map.o vnneutral.o wavefunction.o options.o energy.o grid.o loops.o workmat.o FIRE.o debug.o timing.o

MATH =  cross.o factorial.o cepal.o

MAIN =  fireball.o libFireCore.o denmat.o fermie.o mixer.o anderson2.o sqrtS.o ktransform.o solveH.o

DASSEMBLERS =  Dassemble_2c.o Dassemble_3c.o Dassemble_ca_2c.o Dassemble_ca_3c.o Dassemble_lr.o Dassemble_snxc_on.o Dassemble_olsxc_on.o \
	 Dassemble_olsxc_2c.o Dassemble_olsxc_3c.o Dassemble_snxc_2c.o Dassemble_snxc_3c.o Dassemble_2c_PP.o Dassemble_3c_PP.o Dassemble_ca_olsxc_on.o \
	 Dassemble_ca_snxc_on.o Dassemble_ca_snxc_3c.o Dassemble_ca_olsxc_3c.o Dassemble_ca_snxc_2c.o Dassemble_ca_olsxc_2c.o getforces_mcweda.o \
	 getforces.o

INTERACTIONS =  cl_value.o get_ewald.o smoother.o unocentros.o doscentros.o doscentrosS.o doscentrosPP.o doscentros_vec.o \
	 doscentrosS_vec.o doscentrosPP_vec.o trescentros.o trescentrosS.o Dtrescentros.o DtrescentrosS.o trescentros_vec.o trescentrosS_vec.o \
	 Dtrescentros_vec.o DtrescentrosS_vec.o internalLambda.o tester2c.o

ALLOCATIONS =  allocate_f.o allocate_h.o allocate_neigh.o allocate_rho.o reallocate_f.o reallocate_h.o reallocate_neigh.o \
	 reallocate_rho.o

ASSEMBLERS =  assemble_olsxc_1c.o assemble_2c.o assemble_3c.o assemble_ca_2c.o assemble_3c_PP.o assemble_2c_PP.o assemble_ca_3c.o \
	 assemble_F.o assemble_lr.o assemble_sVNL.o assemble_usr.o buildh.o assemble_olsxc_on.o assemble_olsxc_off.o build_olsxc_on.o \
	 build_olsxc_off.o average_rho.o average_ca_rho.o build_snxc_on.o build_snxc_off.o assemble_snxc_on.o assemble_snxc_off.o \
	 build_ca_snxc_on.o build_ca_olsxc_on.o assemble_h.o assemble_mcweda.o getenergy.o getenergy_mcweda.o assemble_S.o assemble_2c_S.o

INITIALIZERS =  diagnostics.o initcharges.o initconstants.o initboxes.o initneighbors.o make_mu2shell.o make_munu.o make_munuPP.o \
	 initamat.o make_munuS.o get_info_orbital.o initbasics.o

INTERPOLATERS =  buildspline_1d.o interpolate_1d.o interpolate_1d_vec.o interpolate_2d.o interpolate_2d_vec.o recover_2c.o \
	 recover_3c.o recover_PP.o recoverC.o setterp_2d.o recover_S.o buildspline2_1d.o getpsi.o getYlm.o getvna.o

NEIGHBORS =  backnay.o common_neighbors.o find_neigh_max.o mpairnay.o neighbors.o neighbors_pairs.o find_neighPP_max.o neighborsPP.o \
	 common_neighborsPP.o num_neigh_tot.o

READFILES =  append_string.o read_1c.o read_2c.o read_3c.o readdata_2c.o readdata_3c.o readheader_2c.o readheader_3c.o readinfo.o \
	 readparam.o readbasis.o readlvs.o findFdata.o readdata.o readdata_mcweda.o checksum_options.o getsections.o

ROTATIONS =  chooser.o chooserd.o deps2center.o deps3center.o makeDmat.o makeDmatPP.o rotate.o rotated.o rotatedPP.o twister.o \
	 twisterd.o rotatePP.o epsilon.o

GRID =  read_wf.o read_vna.o readgrid.o initgrid.o initgrid_new.o allocate_grid.o project_dens.o project_dens0.o project_orb.o \
	 project_orb_complex.o ew2mesh.o writeout_xsf.o readdata_KS.o initcharges_KS.o initdenmat.o ceperley_alder.o assemble_KS_mat.o \
	 assemble_KS_usr.o assemble_KS_dcc.o assemble_KS_den0.o assemble_KS_den.o assemble_KS_vna.o laplace_fft.o mixer_KS.o

OBJECTS =  $(MODULES)  $(MATH)  $(MAIN)  $(DASSEMBLERS)  $(INTERACTIONS)  $(ALLOCATIONS)  $(ASSEMBLERS)  $(INITIALIZERS)  \
	 $(INTERPOLATERS)  $(NEIGHBORS)  $(READFILES)  $(ROTATIONS)  $(GRID) 

all :
	make fireball.x
	make libFireCore.so

.PHONY : clean veryclean extraclean

clean : 
	rm -f -r core *.o .nfs* rii_files fireball.x.ip*  *.mod ldtmp* *.vo *~ *.il

veryclean :  clean
	rm -f fireball.x libfireball.a libFireCore.so

extraclean : veryclean

fireball.x : $(OBJECTS)
	$(F90) -o fireball.x  $(FFLAGS) $(OBJECTS) $(LFLAGS)

libFireCore.so : $(OBJECTS)
	$(F90) -o libFireCore.so   -shared -fPIC $(FFLAGS) $(OBJECTS) $(LFLAGS)



#**************************************************
#   MODULES
#**************************************************
#====== variant : ''
dimensions.o : ../fortran/MODULES/dimensions.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/MODULES/dimensions.f90
charges.o : ../fortran/MODULES/charges.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/MODULES/charges.f90
configuration.o : ../fortran/MODULES/configuration.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/MODULES/configuration.f90
constants_fireball.o : ../fortran/MODULES/constants_fireball.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/MODULES/constants_fireball.f90
density.o : ../fortran/MODULES/density.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/MODULES/density.f90
forces.o : ../fortran/MODULES/forces.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/MODULES/forces.f90
Fdata3c.o : ../fortran/MODULES/Fdata3c.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/MODULES/Fdata3c.f90
integrals.o : ../fortran/MODULES/integrals.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/MODULES/integrals.f90
interactions.o : ../fortran/MODULES/interactions.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/MODULES/interactions.f90
kpoints.o : ../fortran/MODULES/kpoints.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/MODULES/kpoints.f90
neighbor_map.o : ../fortran/MODULES/neighbor_map.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/MODULES/neighbor_map.f90
vnneutral.o : ../fortran/MODULES/vnneutral.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/MODULES/vnneutral.f90
wavefunction.o : ../fortran/MODULES/wavefunction.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/MODULES/wavefunction.f90
options.o : ../fortran/MODULES/options.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/MODULES/options.f90
energy.o : ../fortran/MODULES/energy.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/MODULES/energy.f90
grid.o : ../fortran/MODULES/grid.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/MODULES/grid.f90
loops.o : ../fortran/MODULES/loops.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/MODULES/loops.f90
workmat.o : ../fortran/MODULES/workmat.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/MODULES/workmat.f90
FIRE.o : ../fortran/MODULES/FIRE.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/MODULES/FIRE.f90
debug.o : ../fortran/MODULES/debug.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/MODULES/debug.f90
timing.o : ../fortran/MODULES/timing.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/MODULES/timing.f90




#**************************************************
#   MATH
#**************************************************
#====== variant : ''
cross.o : ../fortran/MATH/cross.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/MATH/cross.f90
factorial.o : ../fortran/MATH/factorial.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/MATH/factorial.f90
cepal.o : ../fortran/MATH/cepal.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/MATH/cepal.f90




#**************************************************
#   MAIN
#**************************************************
#====== variant : ''
fireball.o : ../fortran/MAIN/fireball.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/MAIN/fireball.f90
libFireCore.o : ../fortran/MAIN/libFireCore.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/MAIN/libFireCore.f90
denmat.o : ../fortran/MAIN/denmat.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/MAIN/denmat.f90
fermie.o : ../fortran/MAIN/fermie.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/MAIN/fermie.f90
mixer.o : ../fortran/MAIN/mixer.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/MAIN/mixer.f90
anderson2.o : ../fortran/MAIN/anderson2.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/MAIN/anderson2.f90
sqrtS.o : ../fortran/MAIN/sqrtS.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/MAIN/sqrtS.f90
ktransform.o : ../fortran/MAIN/ktransform.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/MAIN/ktransform.f90
solveH.o : ../fortran/MAIN/solveH.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/MAIN/solveH.f90




#**************************************************
#   DASSEMBLERS
#**************************************************
#====== variant : ''
Dassemble_2c.o : ../fortran/DASSEMBLERS/Dassemble_2c.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/DASSEMBLERS/Dassemble_2c.f90
Dassemble_3c.o : ../fortran/DASSEMBLERS/Dassemble_3c.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/DASSEMBLERS/Dassemble_3c.f90
Dassemble_ca_2c.o : ../fortran/DASSEMBLERS/Dassemble_ca_2c.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/DASSEMBLERS/Dassemble_ca_2c.f90
Dassemble_ca_3c.o : ../fortran/DASSEMBLERS/Dassemble_ca_3c.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/DASSEMBLERS/Dassemble_ca_3c.f90
Dassemble_lr.o : ../fortran/DASSEMBLERS/Dassemble_lr.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/DASSEMBLERS/Dassemble_lr.f90
Dassemble_snxc_on.o : ../fortran/DASSEMBLERS/Dassemble_snxc_on.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/DASSEMBLERS/Dassemble_snxc_on.f90
Dassemble_olsxc_on.o : ../fortran/DASSEMBLERS/Dassemble_olsxc_on.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/DASSEMBLERS/Dassemble_olsxc_on.f90
Dassemble_olsxc_2c.o : ../fortran/DASSEMBLERS/Dassemble_olsxc_2c.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/DASSEMBLERS/Dassemble_olsxc_2c.f90
Dassemble_olsxc_3c.o : ../fortran/DASSEMBLERS/Dassemble_olsxc_3c.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/DASSEMBLERS/Dassemble_olsxc_3c.f90
Dassemble_snxc_2c.o : ../fortran/DASSEMBLERS/Dassemble_snxc_2c.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/DASSEMBLERS/Dassemble_snxc_2c.f90
Dassemble_snxc_3c.o : ../fortran/DASSEMBLERS/Dassemble_snxc_3c.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/DASSEMBLERS/Dassemble_snxc_3c.f90
Dassemble_2c_PP.o : ../fortran/DASSEMBLERS/Dassemble_2c_PP.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/DASSEMBLERS/Dassemble_2c_PP.f90
Dassemble_3c_PP.o : ../fortran/DASSEMBLERS/Dassemble_3c_PP.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/DASSEMBLERS/Dassemble_3c_PP.f90
Dassemble_ca_olsxc_on.o : ../fortran/DASSEMBLERS/Dassemble_ca_olsxc_on.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/DASSEMBLERS/Dassemble_ca_olsxc_on.f90
Dassemble_ca_snxc_on.o : ../fortran/DASSEMBLERS/Dassemble_ca_snxc_on.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/DASSEMBLERS/Dassemble_ca_snxc_on.f90
Dassemble_ca_snxc_3c.o : ../fortran/DASSEMBLERS/Dassemble_ca_snxc_3c.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/DASSEMBLERS/Dassemble_ca_snxc_3c.f90
Dassemble_ca_olsxc_3c.o : ../fortran/DASSEMBLERS/Dassemble_ca_olsxc_3c.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/DASSEMBLERS/Dassemble_ca_olsxc_3c.f90
Dassemble_ca_snxc_2c.o : ../fortran/DASSEMBLERS/Dassemble_ca_snxc_2c.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/DASSEMBLERS/Dassemble_ca_snxc_2c.f90
Dassemble_ca_olsxc_2c.o : ../fortran/DASSEMBLERS/Dassemble_ca_olsxc_2c.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/DASSEMBLERS/Dassemble_ca_olsxc_2c.f90
getforces_mcweda.o : ../fortran/DASSEMBLERS/getforces_mcweda.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/DASSEMBLERS/getforces_mcweda.f90
getforces.o : ../fortran/DASSEMBLERS/getforces.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/DASSEMBLERS/getforces.f90




#**************************************************
#   INTERACTIONS
#**************************************************
#====== variant : ''
cl_value.o : ../fortran/INTERACTIONS/cl_value.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/INTERACTIONS/cl_value.f90
get_ewald.o : ../fortran/INTERACTIONS/get_ewald.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/INTERACTIONS/get_ewald.f90
smoother.o : ../fortran/INTERACTIONS/smoother.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/INTERACTIONS/smoother.f90
unocentros.o : ../fortran/INTERACTIONS/unocentros.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/INTERACTIONS/unocentros.f90
doscentros.o : ../fortran/INTERACTIONS/doscentros.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/INTERACTIONS/doscentros.f90
doscentrosS.o : ../fortran/INTERACTIONS/doscentrosS.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/INTERACTIONS/doscentrosS.f90
doscentrosPP.o : ../fortran/INTERACTIONS/doscentrosPP.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/INTERACTIONS/doscentrosPP.f90
doscentros_vec.o : ../fortran/INTERACTIONS/doscentros_vec.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/INTERACTIONS/doscentros_vec.f90
doscentrosS_vec.o : ../fortran/INTERACTIONS/doscentrosS_vec.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/INTERACTIONS/doscentrosS_vec.f90
doscentrosPP_vec.o : ../fortran/INTERACTIONS/doscentrosPP_vec.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/INTERACTIONS/doscentrosPP_vec.f90
trescentros.o : ../fortran/INTERACTIONS/trescentros.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/INTERACTIONS/trescentros.f90
trescentrosS.o : ../fortran/INTERACTIONS/trescentrosS.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/INTERACTIONS/trescentrosS.f90
Dtrescentros.o : ../fortran/INTERACTIONS/Dtrescentros.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/INTERACTIONS/Dtrescentros.f90
DtrescentrosS.o : ../fortran/INTERACTIONS/DtrescentrosS.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/INTERACTIONS/DtrescentrosS.f90
trescentros_vec.o : ../fortran/INTERACTIONS/trescentros_vec.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/INTERACTIONS/trescentros_vec.f90
trescentrosS_vec.o : ../fortran/INTERACTIONS/trescentrosS_vec.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/INTERACTIONS/trescentrosS_vec.f90
Dtrescentros_vec.o : ../fortran/INTERACTIONS/Dtrescentros_vec.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/INTERACTIONS/Dtrescentros_vec.f90
DtrescentrosS_vec.o : ../fortran/INTERACTIONS/DtrescentrosS_vec.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/INTERACTIONS/DtrescentrosS_vec.f90
internalLambda.o : ../fortran/INTERACTIONS/internalLambda.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/INTERACTIONS/internalLambda.f90
tester2c.o : ../fortran/INTERACTIONS/tester2c.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/INTERACTIONS/tester2c.f90




#**************************************************
#   ALLOCATIONS
#**************************************************
#====== variant : ''
allocate_f.o : ../fortran/ALLOCATIONS/allocate_f.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/ALLOCATIONS/allocate_f.f90
allocate_h.o : ../fortran/ALLOCATIONS/allocate_h.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/ALLOCATIONS/allocate_h.f90
allocate_neigh.o : ../fortran/ALLOCATIONS/allocate_neigh.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/ALLOCATIONS/allocate_neigh.f90
allocate_rho.o : ../fortran/ALLOCATIONS/allocate_rho.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/ALLOCATIONS/allocate_rho.f90
reallocate_f.o : ../fortran/ALLOCATIONS/reallocate_f.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/ALLOCATIONS/reallocate_f.f90
reallocate_h.o : ../fortran/ALLOCATIONS/reallocate_h.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/ALLOCATIONS/reallocate_h.f90
reallocate_neigh.o : ../fortran/ALLOCATIONS/reallocate_neigh.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/ALLOCATIONS/reallocate_neigh.f90
reallocate_rho.o : ../fortran/ALLOCATIONS/reallocate_rho.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/ALLOCATIONS/reallocate_rho.f90




#**************************************************
#   ASSEMBLERS
#**************************************************
#====== variant : ''
assemble_olsxc_1c.o : ../fortran/ASSEMBLERS/assemble_olsxc_1c.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/ASSEMBLERS/assemble_olsxc_1c.f90
assemble_2c.o : ../fortran/ASSEMBLERS/assemble_2c.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/ASSEMBLERS/assemble_2c.f90
assemble_3c.o : ../fortran/ASSEMBLERS/assemble_3c.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/ASSEMBLERS/assemble_3c.f90
assemble_ca_2c.o : ../fortran/ASSEMBLERS/assemble_ca_2c.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/ASSEMBLERS/assemble_ca_2c.f90
assemble_3c_PP.o : ../fortran/ASSEMBLERS/assemble_3c_PP.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/ASSEMBLERS/assemble_3c_PP.f90
assemble_2c_PP.o : ../fortran/ASSEMBLERS/assemble_2c_PP.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/ASSEMBLERS/assemble_2c_PP.f90
assemble_ca_3c.o : ../fortran/ASSEMBLERS/assemble_ca_3c.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/ASSEMBLERS/assemble_ca_3c.f90
assemble_F.o : ../fortran/ASSEMBLERS/assemble_F.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/ASSEMBLERS/assemble_F.f90
assemble_lr.o : ../fortran/ASSEMBLERS/assemble_lr.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/ASSEMBLERS/assemble_lr.f90
assemble_sVNL.o : ../fortran/ASSEMBLERS/assemble_sVNL.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/ASSEMBLERS/assemble_sVNL.f90
assemble_usr.o : ../fortran/ASSEMBLERS/assemble_usr.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/ASSEMBLERS/assemble_usr.f90
buildh.o : ../fortran/ASSEMBLERS/buildh.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/ASSEMBLERS/buildh.f90
assemble_olsxc_on.o : ../fortran/ASSEMBLERS/assemble_olsxc_on.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/ASSEMBLERS/assemble_olsxc_on.f90
assemble_olsxc_off.o : ../fortran/ASSEMBLERS/assemble_olsxc_off.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/ASSEMBLERS/assemble_olsxc_off.f90
build_olsxc_on.o : ../fortran/ASSEMBLERS/build_olsxc_on.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/ASSEMBLERS/build_olsxc_on.f90
build_olsxc_off.o : ../fortran/ASSEMBLERS/build_olsxc_off.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/ASSEMBLERS/build_olsxc_off.f90
average_rho.o : ../fortran/ASSEMBLERS/average_rho.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/ASSEMBLERS/average_rho.f90
average_ca_rho.o : ../fortran/ASSEMBLERS/average_ca_rho.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/ASSEMBLERS/average_ca_rho.f90
build_snxc_on.o : ../fortran/ASSEMBLERS/build_snxc_on.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/ASSEMBLERS/build_snxc_on.f90
build_snxc_off.o : ../fortran/ASSEMBLERS/build_snxc_off.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/ASSEMBLERS/build_snxc_off.f90
assemble_snxc_on.o : ../fortran/ASSEMBLERS/assemble_snxc_on.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/ASSEMBLERS/assemble_snxc_on.f90
assemble_snxc_off.o : ../fortran/ASSEMBLERS/assemble_snxc_off.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/ASSEMBLERS/assemble_snxc_off.f90
build_ca_snxc_on.o : ../fortran/ASSEMBLERS/build_ca_snxc_on.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/ASSEMBLERS/build_ca_snxc_on.f90
build_ca_olsxc_on.o : ../fortran/ASSEMBLERS/build_ca_olsxc_on.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/ASSEMBLERS/build_ca_olsxc_on.f90
assemble_h.o : ../fortran/ASSEMBLERS/assemble_h.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/ASSEMBLERS/assemble_h.f90
assemble_mcweda.o : ../fortran/ASSEMBLERS/assemble_mcweda.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/ASSEMBLERS/assemble_mcweda.f90
getenergy.o : ../fortran/ASSEMBLERS/getenergy.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/ASSEMBLERS/getenergy.f90
getenergy_mcweda.o : ../fortran/ASSEMBLERS/getenergy_mcweda.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/ASSEMBLERS/getenergy_mcweda.f90
assemble_S.o : ../fortran/ASSEMBLERS/assemble_S.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/ASSEMBLERS/assemble_S.f90
assemble_2c_S.o : ../fortran/ASSEMBLERS/assemble_2c_S.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/ASSEMBLERS/assemble_2c_S.f90




#**************************************************
#   INITIALIZERS
#**************************************************
#====== variant : ''
diagnostics.o : ../fortran/INITIALIZERS/diagnostics.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/INITIALIZERS/diagnostics.f90
initcharges.o : ../fortran/INITIALIZERS/initcharges.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/INITIALIZERS/initcharges.f90
initconstants.o : ../fortran/INITIALIZERS/initconstants.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/INITIALIZERS/initconstants.f90
initboxes.o : ../fortran/INITIALIZERS/initboxes.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/INITIALIZERS/initboxes.f90
initneighbors.o : ../fortran/INITIALIZERS/initneighbors.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/INITIALIZERS/initneighbors.f90
make_mu2shell.o : ../fortran/INITIALIZERS/make_mu2shell.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/INITIALIZERS/make_mu2shell.f90
make_munu.o : ../fortran/INITIALIZERS/make_munu.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/INITIALIZERS/make_munu.f90
make_munuPP.o : ../fortran/INITIALIZERS/make_munuPP.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/INITIALIZERS/make_munuPP.f90
initamat.o : ../fortran/INITIALIZERS/initamat.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/INITIALIZERS/initamat.f90
make_munuS.o : ../fortran/INITIALIZERS/make_munuS.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/INITIALIZERS/make_munuS.f90
get_info_orbital.o : ../fortran/INITIALIZERS/get_info_orbital.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/INITIALIZERS/get_info_orbital.f90
initbasics.o : ../fortran/INITIALIZERS/initbasics.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/INITIALIZERS/initbasics.f90




#**************************************************
#   INTERPOLATERS
#**************************************************
#====== variant : ''
buildspline_1d.o : ../fortran/INTERPOLATERS/buildspline_1d.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/INTERPOLATERS/buildspline_1d.f90
interpolate_1d.o : ../fortran/INTERPOLATERS/interpolate_1d.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/INTERPOLATERS/interpolate_1d.f90
interpolate_1d_vec.o : ../fortran/INTERPOLATERS/interpolate_1d_vec.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/INTERPOLATERS/interpolate_1d_vec.f90
interpolate_2d.o : ../fortran/INTERPOLATERS/interpolate_2d.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/INTERPOLATERS/interpolate_2d.f90
interpolate_2d_vec.o : ../fortran/INTERPOLATERS/interpolate_2d_vec.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/INTERPOLATERS/interpolate_2d_vec.f90
recover_2c.o : ../fortran/INTERPOLATERS/recover_2c.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/INTERPOLATERS/recover_2c.f90
recover_3c.o : ../fortran/INTERPOLATERS/recover_3c.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/INTERPOLATERS/recover_3c.f90
recover_PP.o : ../fortran/INTERPOLATERS/recover_PP.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/INTERPOLATERS/recover_PP.f90
recoverC.o : ../fortran/INTERPOLATERS/recoverC.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/INTERPOLATERS/recoverC.f90
setterp_2d.o : ../fortran/INTERPOLATERS/setterp_2d.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/INTERPOLATERS/setterp_2d.f90
recover_S.o : ../fortran/INTERPOLATERS/recover_S.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/INTERPOLATERS/recover_S.f90
buildspline2_1d.o : ../fortran/INTERPOLATERS/buildspline2_1d.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/INTERPOLATERS/buildspline2_1d.f90
getpsi.o : ../fortran/INTERPOLATERS/getpsi.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/INTERPOLATERS/getpsi.f90
getYlm.o : ../fortran/INTERPOLATERS/getYlm.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/INTERPOLATERS/getYlm.f90
getvna.o : ../fortran/INTERPOLATERS/getvna.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/INTERPOLATERS/getvna.f90




#**************************************************
#   NEIGHBORS
#**************************************************
#====== variant : ''
backnay.o : ../fortran/NEIGHBORS/backnay.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/NEIGHBORS/backnay.f90
common_neighbors.o : ../fortran/NEIGHBORS/common_neighbors.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/NEIGHBORS/common_neighbors.f90
find_neigh_max.o : ../fortran/NEIGHBORS/find_neigh_max.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/NEIGHBORS/find_neigh_max.f90
mpairnay.o : ../fortran/NEIGHBORS/mpairnay.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/NEIGHBORS/mpairnay.f90
neighbors.o : ../fortran/NEIGHBORS/neighbors.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/NEIGHBORS/neighbors.f90
neighbors_pairs.o : ../fortran/NEIGHBORS/neighbors_pairs.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/NEIGHBORS/neighbors_pairs.f90
find_neighPP_max.o : ../fortran/NEIGHBORS/find_neighPP_max.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/NEIGHBORS/find_neighPP_max.f90
neighborsPP.o : ../fortran/NEIGHBORS/neighborsPP.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/NEIGHBORS/neighborsPP.f90
common_neighborsPP.o : ../fortran/NEIGHBORS/common_neighborsPP.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/NEIGHBORS/common_neighborsPP.f90
num_neigh_tot.o : ../fortran/NEIGHBORS/num_neigh_tot.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/NEIGHBORS/num_neigh_tot.f90




#**************************************************
#   READFILES
#**************************************************
#====== variant : ''
append_string.o : ../fortran/READFILES/append_string.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/READFILES/append_string.f90
read_1c.o : ../fortran/READFILES/read_1c.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/READFILES/read_1c.f90
read_2c.o : ../fortran/READFILES/read_2c.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/READFILES/read_2c.f90
read_3c.o : ../fortran/READFILES/read_3c.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/READFILES/read_3c.f90
readdata_2c.o : ../fortran/READFILES/readdata_2c.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/READFILES/readdata_2c.f90
readdata_3c.o : ../fortran/READFILES/readdata_3c.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/READFILES/readdata_3c.f90
readheader_2c.o : ../fortran/READFILES/readheader_2c.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/READFILES/readheader_2c.f90
readheader_3c.o : ../fortran/READFILES/readheader_3c.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/READFILES/readheader_3c.f90
readinfo.o : ../fortran/READFILES/readinfo.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/READFILES/readinfo.f90
readparam.o : ../fortran/READFILES/readparam.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/READFILES/readparam.f90
readbasis.o : ../fortran/READFILES/readbasis.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/READFILES/readbasis.f90
readlvs.o : ../fortran/READFILES/readlvs.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/READFILES/readlvs.f90
findFdata.o : ../fortran/READFILES/findFdata.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/READFILES/findFdata.f90
readdata.o : ../fortran/READFILES/readdata.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/READFILES/readdata.f90
readdata_mcweda.o : ../fortran/READFILES/readdata_mcweda.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/READFILES/readdata_mcweda.f90
checksum_options.o : ../fortran/READFILES/checksum_options.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/READFILES/checksum_options.f90
getsections.o : ../fortran/READFILES/getsections.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/READFILES/getsections.f90




#**************************************************
#   ROTATIONS
#**************************************************
#====== variant : ''
chooser.o : ../fortran/ROTATIONS/chooser.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/ROTATIONS/chooser.f90
chooserd.o : ../fortran/ROTATIONS/chooserd.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/ROTATIONS/chooserd.f90
deps2center.o : ../fortran/ROTATIONS/deps2center.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/ROTATIONS/deps2center.f90
deps3center.o : ../fortran/ROTATIONS/deps3center.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/ROTATIONS/deps3center.f90
makeDmat.o : ../fortran/ROTATIONS/makeDmat.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/ROTATIONS/makeDmat.f90
makeDmatPP.o : ../fortran/ROTATIONS/makeDmatPP.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/ROTATIONS/makeDmatPP.f90
rotate.o : ../fortran/ROTATIONS/rotate.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/ROTATIONS/rotate.f90
rotated.o : ../fortran/ROTATIONS/rotated.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/ROTATIONS/rotated.f90
rotatedPP.o : ../fortran/ROTATIONS/rotatedPP.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/ROTATIONS/rotatedPP.f90
twister.o : ../fortran/ROTATIONS/twister.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/ROTATIONS/twister.f90
twisterd.o : ../fortran/ROTATIONS/twisterd.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/ROTATIONS/twisterd.f90
rotatePP.o : ../fortran/ROTATIONS/rotatePP.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/ROTATIONS/rotatePP.f90
epsilon.o : ../fortran/ROTATIONS/epsilon.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/ROTATIONS/epsilon.f90




#**************************************************
#   GRID
#**************************************************
#====== variant : ''
read_wf.o : ../fortran/GRID/read_wf.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/GRID/read_wf.f90
read_vna.o : ../fortran/GRID/read_vna.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/GRID/read_vna.f90
readgrid.o : ../fortran/GRID/readgrid.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/GRID/readgrid.f90
initgrid.o : ../fortran/GRID/initgrid.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/GRID/initgrid.f90
initgrid_new.o : ../fortran/GRID/initgrid_new.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/GRID/initgrid_new.f90
allocate_grid.o : ../fortran/GRID/allocate_grid.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/GRID/allocate_grid.f90
project_dens.o : ../fortran/GRID/project_dens.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/GRID/project_dens.f90
project_dens0.o : ../fortran/GRID/project_dens0.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/GRID/project_dens0.f90
project_orb.o : ../fortran/GRID/project_orb.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/GRID/project_orb.f90
project_orb_complex.o : ../fortran/GRID/project_orb_complex.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/GRID/project_orb_complex.f90
ew2mesh.o : ../fortran/GRID/ew2mesh.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/GRID/ew2mesh.f90
writeout_xsf.o : ../fortran/GRID/writeout_xsf.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/GRID/writeout_xsf.f90
readdata_KS.o : ../fortran/GRID/readdata_KS.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/GRID/readdata_KS.f90
initcharges_KS.o : ../fortran/GRID/initcharges_KS.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/GRID/initcharges_KS.f90
initdenmat.o : ../fortran/GRID/initdenmat.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/GRID/initdenmat.f90
ceperley_alder.o : ../fortran/GRID/ceperley_alder.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/GRID/ceperley_alder.f90
assemble_KS_mat.o : ../fortran/GRID/assemble_KS_mat.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/GRID/assemble_KS_mat.f90
assemble_KS_usr.o : ../fortran/GRID/assemble_KS_usr.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/GRID/assemble_KS_usr.f90
assemble_KS_dcc.o : ../fortran/GRID/assemble_KS_dcc.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/GRID/assemble_KS_dcc.f90
assemble_KS_den0.o : ../fortran/GRID/assemble_KS_den0.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/GRID/assemble_KS_den0.f90
assemble_KS_den.o : ../fortran/GRID/assemble_KS_den.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/GRID/assemble_KS_den.f90
assemble_KS_vna.o : ../fortran/GRID/assemble_KS_vna.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/GRID/assemble_KS_vna.f90
laplace_fft.o : ../fortran/GRID/laplace_fft.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/GRID/laplace_fft.f90
mixer_KS.o : ../fortran/GRID/mixer_KS.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c ../fortran/GRID/mixer_KS.f90


