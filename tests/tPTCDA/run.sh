#!/bin/bash

root_firecore=/home/niko/work/GRIDFF/FireCore                          # put here your actual path to main FireCore folder
lammps_exe=/home/niko/work/src/lammps/lammps-29Aug2024/src/lmp_smiggol # put here your actual path to LAMMPS executable

mol=PTCDA

# generate confs
nconfs=41
runs="distort xscan yscan zscan"
dist_max=0.05 # Ang
xscanmin=-2.0 # Ang, xscanmax=2.0 # Ang
yscanmin=$xscanmin # yscanmax=$xscanmax
zscanmin=-0.3 # Ang, zscanmax=3.7 # Ang
dscan=0.1 # Ang

# test interactions
test_key1s=(
    "dist" #"01-bonds"			       
    "dist" #"02-angles"			       
    "dist" #"03-dihedrals"		       
    "dist" #"04-inversions"		       
    "dist" #"05-lj" 			       
    "dist" #"06-bonds+lj_excl_1-2_sub_noclamp" 
    "dist" #"07-bonds+lj_excl_1-2_sub_clamp"   
    "dist" #"08-bonds+lj_excl_1-2_ng4"	       
    "dist" #"09-intra_excl_sub_noclamp"	       
    "dist" #"10-intra_excl_sub_clamp"	       
    "dist" #"11-intra_excl_ng4_noclamp"	       
    "dist" #"12-intra_excl_ng4_clamp"	       
    "both" #"13-sub"			       
    "nacl" #"14-all_excl_sub_noclamp"	       
    "nacl" #"15-all_excl_sub_clamp"	       
    "nacl" #"16-all_excl_ng4_noclamp"	       
    "nacl" #"17-all_excl_ng4_clamp"            
)
test_key2s=(
    "01-bonds"
    "02-angles"
    "03-dihedrals"
    "04-inversions"
    "05-lj" 
    "06-bonds+lj_excl_1-2_sub_noclamp"
    "07-bonds+lj_excl_1-2_sub_clamp"
    "08-bonds+lj_excl_1-2_ng4"
    "09-intra_excl_sub_noclamp"
    "10-intra_excl_sub_clamp"
    "11-intra_excl_ng4_noclamp"
    "12-intra_excl_ng4_clamp"
    "13-sub"
    "14-all_excl_sub_noclamp"
    "15-all_excl_sub_clamp"
    "16-all_excl_ng4_noclamp"
    "17-all_excl_ng4_clamp"            
)
test_key3s=(
    "Bonds only"                                                                       #"01-bonds"			       
    "Angles only"									#"02-angles"			       
    "Dihedrals only"									#"03-dihedrals"		       
    "Inversions only"									#"04-inversions"		       
    "Intramolecular LJ only"								#"05-lj" 			       
    "Bonds+LJ excluding 1-2 interactions with subtraction"				#"06-bonds+lj_excl_1-2_sub_noclamp" 
    "Bonds+LJ excluding 1-2 interactions with subtraction/clamping"			#"07-bonds+lj_excl_1-2_sub_clamp"   
    "Bonds+LJ excluding 1-2 interactions with ng4"					#"08-bonds+lj_excl_1-2_ng4"	       
    "All intramolecular (exclusion of 1-2 and 1-3 with subtraction)"			#"09-intra_excl_sub_noclamp"	       
    "All intramolecular (exclusion of 1-2 and 1-3 with subtraction/clamping)"		#"10-intra_excl_sub_clamp"	       
    "All intramolecular (exclusion of 1-2 with ng4 and 1-3 with subtraction)"		#"11-intra_excl_ng4_noclamp"	       
    "All intramolecular (exclusion of 1-2 with ng4 and 1-3 with subtraction/clamping)"	#"12-intra_excl_ng4_clamp"	       
    "All interaction with the substrate only"						#"13-sub"			       
    "All interactions (exclusion of 1-2 and 1-3 with subtraction)"			#"14-all_excl_sub_noclamp"	       
    "All interactions (exclusion of 1-2 and 1-3 with subtraction/clamping)"		#"15-all_excl_sub_clamp"	       
    "All interactions (exclusion of 1-2 with ng4 and 1-3 with subtraction)"		#"16-all_excl_ng4_noclamp"	       
    "All interactions (exclusion of 1-2 with ng4 and 1-3 with subtraction/clamping)"	#"17-all_excl_ng4_clamp"            
)
ntests=${#test_key1s[@]}

##########################################################

# generate confs
mkdir -p confs
cd confs
gfortran -Wall -Wextra -Wconversion -pedantic -O -fcheck=all -g -fbacktrace -ffpe-trap=zero,overflow,underflow,invalid ../f90/shift.f90 -o shift.x
ln -sf ../inputs/$mol.data
for run in $runs ; do
    echo generate confs $run
    if [ $run == 'distort' ] ; then shift=$dist_max
    elif [ $run == 'xscan' ] ; then shift=`echo "$xscanmin - $dscan" | bc -l`
    elif [ $run == 'yscan' ] ; then shift=`echo "$yscanmin - $dscan" | bc -l`
    elif [ $run == 'zscan' ] ; then shift=`echo "$zscanmin - $dscan" | bc -l`
    fi
    for (( ic=1 ; ic <= nconfs ; ic++ )) ; do
	if [ $run != 'distort' ] ; then shift=`echo "$shift + $dscan" | bc -l` ; fi
	echo "$mol.data $mol.$run$ic.data $mol.$run$ic.xyz sub.xyz $run $shift" | ./shift.x
    done
done
rm shift.x $mol.data
cd ..

# run LAMMPS for reference
mkdir -p ref
cd ref
for (( it=0 ; it < ntests ; it++ )) ; do
    t1=${test_key1s[$it]}
    t2=${test_key2s[$it]}
    mkdir -p $t2
    cd $t2
    for run in $runs ; do
	if [ $t1 == 'dist' ] && [ ${run:1:4} == 'scan' ] ; then continue ; fi
	echo reference $t2 $run
	for (( ic=1 ; ic <= nconfs ; ic++ )) ; do
	    mkdir -p $run$ic
	    cd $run$ic
	    ln -sf ../../../confs/$mol.$run$ic.data init.data
	    ln -sf ../../../inputs/$t2.in lammps.in
	    $lammps_exe -in lammps.in -log lammps.log 1> /dev/null
	    natoms=`grep atoms init.data | awk '{print $1}'`
	    tail -$natoms traj.lammpstrj > f_lammps.txt
	    grep -B1 Loop lammps.log | head -1 > e_lammps.txt
	    cd ..
	done
    done
    cd ..
done
cd ..

# run FireCore
cd $root_firecore/cpp/Build/libs/Molecular
make -j MMFF_lib 1> /dev/null
cd - 1> /dev/null
mkdir -p firecore
cd firecore
for (( it=0 ; it < ntests ; it++ )) ; do
    t1=${test_key1s[$it]}
    t2=${test_key2s[$it]}
    t3=${test_key3s[$it]}
    mkdir -p $t2
    cd $t2
    for run in $runs ; do
	if [ $t1 == 'dist' ] && [ ${run:1:4} == 'scan' ] ; then continue ; fi
	for (( ic=1 ; ic <= nconfs ; ic++ )) ; do
	    echo firecore $t2 $run conf$ic
	    mkdir -p $run$ic
	    cd $run$ic
	    ln -sf $root_firecore/cpp/common_resources data
	    ln -sf $root_firecore/cpp/common_resources common_resources 
	    ln -sf ../../../confs/$mol.$run$ic.xyz mol.xyz
	    ln -sf ../../../confs/sub.xyz
	    x=`python3 ../../../run.py ${t2:0:2} &> firecore.log`
	    cd ..
	done
    done
    cd ..
done
cd ..

# compare & plot
cd firecore
gfortran -Wall -Wextra -Wconversion -pedantic -O -fcheck=all -g -fbacktrace -ffpe-trap=zero,overflow,underflow,invalid ../f90/compare.f90 -o compare.x
gfortran -Wall -Wextra -Wconversion -pedantic -O -fcheck=all -g -fbacktrace -ffpe-trap=zero,overflow,underflow,invalid ../f90/subtract.f90 -o subtract.x
for (( it=0 ; it < ntests ; it++ )) ; do
    t1=${test_key1s[$it]}
    t2=${test_key2s[$it]}
    t3=${test_key3s[$it]}
    mkdir -p $t2
    cd $t2
    for run in $runs ; do
	if [ $t1 == 'dist' ] && [ ${run:1:4} == 'scan' ] ; then continue ; fi
	echo plot $t2 $run
	rm -f ${run}_DIFF_FORCES ${run}_DIFF_ENERGY
	for (( ic=1 ; ic <= nconfs ; ic++ )) ; do
	    mkdir -p $run$ic
	    cd $run$ic
	    rm -f ?_lammps.txt
	    if [ $t1 == 'nacl' ] || [ $t1 == 'both' ] ; then
		ln -sf ../../../ref/$t2/$run$ic/f_lammps.txt f_current.txt
		ln -sf ../../../ref/$t2/$run$ic/e_lammps.txt e_current.txt
		ln -sf ../../../confs/$mol.$run$ic.data init.data
		natoms=`grep atoms init.data | awk '{print $1}'`
		ln -sf ../../../inputs/coul_nacl.in lammps.in
		$lammps_exe -in lammps.in -log lammps.log 1> /dev/null
		tail -$natoms traj.lammpstrj > f_offset.txt
		grep -B1 Loop lammps.log | head -1 > e_offset.txt
		echo $natoms | ../../subtract.x
		rm ?_current.txt init.data lammps.in lammps.log traj.lammpstrj ?_offset.txt
		if [ $t1 == 'both' ] ; then
		    mv f_lammps.txt f_current.txt
		    mv e_lammps.txt e_current.txt
		    ln -sf ../../../confs/$mol.$run$ic.data init.data
		    ln -sf ../../../inputs/coul_ptcda.in lammps.in
		    $lammps_exe -in lammps.in -log lammps.log 1> /dev/null
		    tail -$natoms traj.lammpstrj > f_offset.txt
		    grep -B1 Loop lammps.log | head -1 > e_offset.txt
		    echo $natoms | ../../subtract.x
		    rm ?_current.txt init.data lammps.in lammps.log traj.lammpstrj ?_offset.txt
		fi
	    else
		ln -sf ../../../ref/$t2/$run$ic/?_lammps.txt .
	    fi
	    natoms=`head -1 mol.xyz | awk '{print $1}'`
	    echo $natoms | ../../compare.x > DIFF
	    grep FORCE DIFF >> ../${run}_DIFF_FORCES
	    grep ENERGY DIFF >> ../${run}_DIFF_ENERGY
	    cd ..
	done
	echo "set terminal png" > x.gp
	echo "set output '../../$t2.diff.$run.png'" >> x.gp
	echo "set title '$t3'" >> x.gp
	echo "set key bottom right" >> x.gp
	echo "set xrange [0:$((natoms*3*nconfs-1))]" >> x.gp
	echo "unset xtics" >> x.gp
	echo "set ylabel 'Force difference (eV/Ang)'" >> x.gp
	echo "set ytics nomirror" >> x.gp
	echo "set y2label 'Energy difference (eV)'" >> x.gp
	echo "set y2tics nomirror" >> x.gp
	echo "set logscale y" >> x.gp
	echo "set logscale y2" >> x.gp
	echo "plot '${run}_DIFF_FORCES' u 0:6 w p ps 1 title 'forces', "'\' >> x.gp
	echo "     '${run}_DIFF_ENERGY' u ("'$0'"*3*"$natoms"):4 w p pt 5 axes x1y2 title 'energies'" >> x.gp
	gnuplot x.gp
	if [ $run != 'distort' ] ; then
	    echo "set terminal png" > x.gp
	    echo "set output '../../$t2.energy.$run.png'" >> x.gp
	    echo "set title '$t3'" >> x.gp
	    echo "set key bottom right" >> x.gp
	    echo "set xlabel 'Shift $run (Ang)'" >> x.gp
	    echo "set ylabel 'Energy (eV)'" >> x.gp
	    echo "plot '${run}_DIFF_ENERGY' u ("'$0'"*"$dscan"):2 w lp title 'LAMMPS', "'\' >> x.gp
	    echo "     '${run}_DIFF_ENERGY' u ("'$0'"*"$dscan"):3 w l title 'FireCore'" >> x.gp
	    gnuplot x.gp
	fi
    done
    rm x.gp
    cd ..
done
rm compare.x subtract.x
cd ..

exit
