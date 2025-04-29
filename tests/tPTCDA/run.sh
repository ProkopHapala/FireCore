#!/bin/bash

root_firecore=/home/indranil/git/FireCore  
##lammps_exe=/home/niko/work/src/lammps/lammps-29Aug2024/src/lmp_smiggol
lammps_exe=/home/indranil/src/lammps-29Aug2024/src/lmp_mpi
##########################################################

mol=PTCDA

# confs
nconfs=41
runs="distort xscan yscan zscan"
dist_max=0.05 # Ang
xscanmin=-2.0 # Ang, xscanmax=2.0 # Ang
yscanmin=$xscanmin # yscanmax=$xscanmax
zscanmin=-1.0 # Ang, zscanmax=3.0 # Ang
dscan=0.1 # Ang

# interactions
test_key2s=("01-bonds"   "02-angles"   "03-dihedrals"   "04-inversions"   "13-all_ng4_noclamp" "20-nb_sub_only"             "21-all")
test_key3s=("Bonds only" "Angles only" "Dihedrals only" "Inversions only" "All intramolecular" "Interaction with substrate" "Everything")
ntests=${#test_key2s[@]}

# flags
confs_flag=true
#confs_flag=false
ref_flag=true
#ref_flag=false
run_flag=true
#run_flag=false

##########################################################

# generate confs
if $confs_flag ; then
    rm -rf confs
    mkdir -p confs
    cd confs
    gfortran -Wall -Wextra -Wconversion -pedantic -O -fcheck=all -g -fbacktrace -ffpe-trap=zero,overflow,underflow,invalid ../shift.f90 -o shift.x
    ln -sf ../$mol.data
    for run in $runs ; do
	echo generate confs $run
	if [ "$run" == 'xscan' ] ; then
	    shift=`echo "$xscanmin - $dscan" | bc -l`
	elif [ "$run" == 'yscan' ] ; then
	    shift=`echo "$yscanmin - $dscan" | bc -l`
	elif [ "$run" == 'zscan' ] ; then
	    shift=`echo "$zscanmin - $dscan" | bc -l`
	fi
	for (( ic=1 ; ic <= nconfs ; ic++ )) ; do
	    if [ "$run" == 'distort' ] ; then
		echo "$mol.data $mol.$run$ic.data $mol.$run$ic.xyz sub.xyz $run $dist_max" | ./shift.x
	    else
		shift=`echo "$shift + $dscan" | bc -l`
		echo "$mol.data $mol.$run$ic.data $mol.$run$ic.xyz sub.xyz $run $shift" | ./shift.x
	    fi
	done
    done
    rm shift.x $mol.data
    cd ..
fi

# run LAMMPS for reference
if $ref_flag ; then
    rm -rf ref
    mkdir -p ref
    cd ref
    for (( it=0 ; it < ntests ; it++ )) ; do
	t2=${test_key2s[$it]}
	mkdir -p $t2
	cd $t2
        for run in $runs ; do
	    if ( [ ${t2:0:2} == '01' ] || [ ${t2:0:2} == '02' ] || [ ${t2:0:2} == '03' ] || [ ${t2:0:2} == '04' ] || [ ${t2:0:2} == '13' ] ) && [ ${run:1:4} == 'scan' ] ; then continue ; fi
	    echo reference $t2 $run
	    for (( ic=1 ; ic <= nconfs ; ic++ )) ; do
		mkdir -p $run$ic
		cd $run$ic
		ln -sf ../../../confs/$mol.$run$ic.data init.data
		ln -sf ../../../$t2.in lammps.in
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
fi

# run FireCore
if $run_flag ; then
    cd $root_firecore/cpp/Build/libs/Molecular
    make -j MMFF_lib 1> /dev/null
    cd - 1> /dev/null
    rm -rf firecore *.png
    mkdir -p firecore
    cd firecore
    gfortran -Wall -Wextra -Wconversion -pedantic -O -fcheck=all -g -fbacktrace -ffpe-trap=zero,overflow,underflow,invalid ../compare.f90 -o compare.x
    for (( it=0 ; it < ntests ; it++ )) ; do
	t2=${test_key2s[$it]}
	t3=${test_key3s[$it]}
	mkdir -p $t2
	cd $t2
	for run in $runs ; do
	    if ( [ ${t2:0:2} == '01' ] || [ ${t2:0:2} == '02' ] || [ ${t2:0:2} == '03' ] || [ ${t2:0:2} == '04' ] || [ ${t2:0:2} == '13' ] ) && [ ${run:1:4} == 'scan' ] ; then continue ; fi
	    rm -f ${run}_DIFF_FORCES ${run}_DIFF_ENERGY
	    for (( ic=1 ; ic <= nconfs ; ic++ )) ; do
		echo firecore $t2 $run conf$ic
		mkdir -p $run$ic
		cd $run$ic
		ln -sf ../../../ref/$t2/$run$ic/?_lammps.txt .
		ln -sf $root_firecore/cpp/common_resources data
		ln -sf $root_firecore/cpp/common_resources common_resources 
		ln -sf ../../../confs/$mol.$run$ic.xyz mol.xyz
		ln -sf ../../../confs/sub.xyz
		python3 ../../../firecore.py ${t2:0:2} | tee firecore.log #&> firecore.log
		natoms=`head -1 mol.xyz | awk '{print $1}'`
		echo $natoms | ../../compare.x > DIFF
		grep FORCE DIFF >> ../${run}_DIFF_FORCES
		grep ENERGY DIFF >> ../${run}_DIFF_ENERGY
		cd ..
	    done
	    echo "set terminal png" > x.gp
	    echo "set output '../../$t2.diff.$run.png'" >> x.gp
	    echo "set title '$t3'" >> x.gp
	    echo "set xrange [0:$((natoms*3*nconfs-1))]" >> x.gp
	    echo "unset xtics" >> x.gp
	    echo "set ylabel 'Force difference (eV/Ang)'" >> x.gp
	    echo "set ytics nomirror" >> x.gp
	    echo "set y2label 'Energy difference (eV)'" >> x.gp
	    echo "set y2tics nomirror" >> x.gp
	    echo "set logscale y" >> x.gp
	    echo "set logscale y2" >> x.gp
	    echo "plot '${run}_DIFF_FORCES' u 0:6 w p ps 1 title 'forces', "'\' >> x.gp
	    echo "'${run}_DIFF_ENERGY' u ("'$0'"*3*"$natoms"):4 w p pt 5 axes x1y2 title 'energies'" >> x.gp
	    gnuplot x.gp
	    if [ "$run" != 'distort' ] ; then
		echo "set terminal png" > x.gp
		echo "set output '../../$t2.energy.$run.png'" >> x.gp
		echo "set key bottom right" >> x.gp
		echo "set title '$t3'" >> x.gp
		echo "set xlabel 'Shift (Ang)'" >> x.gp
		echo "set ylabel 'Energy (eV)'" >> x.gp
		echo "plot '${run}_DIFF_ENERGY' u ("'$0'"*0.1):2 w lp title 'LAMMPS', "'\' >> x.gp
		echo "     '${run}_DIFF_ENERGY' u ("'$0'"*0.1):3 w l title 'FireCore'" >> x.gp
		gnuplot x.gp
	    fi
	done
	rm x.gp
	cd ..
    done
    rm compare.x
    cd ..
fi

exit
