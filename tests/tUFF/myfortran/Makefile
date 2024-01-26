FC = gfortran

DEBUG ?= 1
ifeq ($(DEBUG), 1)
    FFLAGS = -I. -ffree-line-length-none -Wall -Wextra -Wconversion -pedantic -O -fcheck=all -g -fbacktrace -ffpe-trap=zero,overflow,underflow,invalid
else
    FFLAGS = -I. -ffree-line-length-none -O3
endif

EXE = myfortran.x

OBJECTS = modules.o \
	main.o \
	read_uff_params.o \
	read_conf.o cross_product.o \
	find_neigh.o inverse.o calc_dist.o calc_dist_eq.o \
	find_bond.o trivial_cases.o special_case_nitro.o special_case_amide.o special_case_cumulene.o \
	tree_walk.o check_all.o check_sat.o check_done.o check_atom.o \
	find_rings.o check_aroma.o set_aroma.o \
	simple_rule.o \
	assign_O2.o assign_O3.o assign_N2.o assign_N3.o assign_C2.o assign_C1.o \
	fix_saturation.o \
	conjugation.o \
	checks.o \
	assign_vdw.o \
	assign_bonds.o find_I.o calc_r_IJ.o \
	assign_angles.o \
	assign_torsions.o \
	assign_inversions.o \
	write_data.o 

all: $(OBJECTS)
	$(FC) -o $(EXE) $(FFLAGS) $(OBJECTS)

clean:
	-rm -f *__genmod.f90 *.o *.mod $(EXE)

%.o: %.f90
	$(FC) $(FFLAGS) -c $<
