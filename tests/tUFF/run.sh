#!/bin/bash

mol=PTCDA
natoms=38
nconf=20

#t1s=("01" "02" "03" "04" "05" "06" "07" "08" "09" "10" "11" "12" "13" "14")
#t2s=("01-bonds" "02-angles" "03-dihedrals" "04-inversions" "05-nb" \
#		"06-nb_excluding_bonds_subtract_noclamp" "07-nb_excluding_bonds_subtract_clamp" \
#                "08-nb_excluding_bonds_ng4" "09-nb_excluding_angles_subtract_noclamp" \
#		"10-nb_excluding_angles_subtract_clamp" "11-all_subtract_noclamp" \
#		"12-all_subtract_clamp" "13-all_ng4_noclamp" "14-all_ng4_clamp")
#t3s=("Bonds only" "Angles only" "Dihedrals only" "Inversions only" "Non-bonded only" \
#		  "Non-bonded only + exclusion of 1-2 with subtraction" "Non-bonded only + exclusion of 1-2 with subtraction and clamping" \
#		  "Non-bonded only + exclusion of 1-2 with ng4" "Non-bonded only + exclusion of 1-3 with subtraction" \
#		  "Non-bonded only + exclusion of 1-3 with subtraction and clamping" "UFF + exclusion with subtraction" \
#		  "UFF + exclusion with subtraction and clamping" "UFF + exclusion with ng4" "UFF + exclusion with ng4 and clamping")
t1s=("01" "02" "03" "04" "05" "06" "07" "08" "09" "11" "12" "13" "14")
t2s=("01-bonds" "02-angles" "03-dihedrals" "04-inversions" "05-nb" \
		"06-nb_excluding_bonds_subtract_noclamp" "07-nb_excluding_bonds_subtract_clamp" \
                "08-nb_excluding_bonds_ng4" "09-nb_excluding_angles_subtract_noclamp" \
		"11-all_subtract_noclamp" \
		"12-all_subtract_clamp" "13-all_ng4_noclamp" "14-all_ng4_clamp")
t3s=("Bonds only" "Angles only" "Dihedrals only" "Inversions only" "Non-bonded only" \
		  "Non-bonded only + exclusion of 1-2 with subtraction" "Non-bonded only + exclusion of 1-2 with subtraction and clamping" \
		  "Non-bonded only + exclusion of 1-2 with ng4" "Non-bonded only + exclusion of 1-3 with subtraction" \
		  "UFF + exclusion with subtraction" \
		  "UFF + exclusion with subtraction and clamping" "UFF + exclusion with ng4" "UFF + exclusion with ng4 and clamping")
nt=${#t1s[@]}

cd ../../cpp/Build/libs/Molecular
make -j4 MMFF_lib 1> /dev/null
cd - 1> /dev/null

gfortran -Wall compare.f90 -o compare.x

for (( it=0 ; it < nt ; it++ )) ; do
    t1=${t1s[$it]}
    t2=${t2s[$it]}
    t3=${t3s[$it]}
    echo test $t2
    mkdir -p $t2
    cd $t2
    rm -f DIFF_FORCES DIFF_ENERGY
    for (( ic=1 ; ic <= nconf ; ic++ )) ; do
	echo "  conf$ic"
	rm -rf conf$ic
	mkdir -p conf$ic
	cd conf$ic
	ln -sf ../../ref_lammps/$t2/conf$ic/?_lammps.txt .
	ln -sf ../../../../cpp/common_resources data
	ln -sf ../../../../cpp/common_resources common_resources 
	ln -sf ../../ref_confs/PTCDA.$ic.xyz PTCDA.xyz
	python3 ../../run.py $t1 &> firecore.log
	echo $natoms | ../../compare.x > DIFF
	grep FORCE DIFF >> ../DIFF_FORCES
	grep ENERGY DIFF >> ../DIFF_ENERGY
	cd ..
    done
    cd ..
    echo "set terminal png" > x.gp
    echo "set output '$t2.png'" >> x.gp
    echo "set title '$t3'" >> x.gp
    echo "set key bottom right" >> x.gp
    echo "set xlabel 'Entry'" >> x.gp
    echo "set xrange [0:$((natoms*3*nconf-1))]" >> x.gp
    echo "set ylabel 'Relative difference'" >> x.gp
    echo "set yrange [1e-16:*]" >> x.gp
    echo "set logscale y" >> x.gp
    echo "plot '$t2/DIFF_FORCES' u 0:6 w p ps 1 title 'forces', '$t2/DIFF_ENERGY' u ("'$0'"*3*"$natoms"):4 w p pt 5 title 'energies'" >> x.gp
    gnuplot x.gp
done

rm compare.x x.gp

exit
