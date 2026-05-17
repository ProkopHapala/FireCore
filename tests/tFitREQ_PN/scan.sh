#!/bin/bash

### inputs
files_inputs='H2O-A1_H2O-D1-y.xyz'
#files_inputs='H2O-A1_H2O-D1-y.xyz H2O-A1_H2O-D1-z.xyz H2O-D1_H2O-A1-y.xyz H2O-D1_H2O-A1-z.xyz'

### kMorse
#kMorses='1.6 1.7 1.8 -1' # value of alpha for Morse and Buck interactions, a negative value means that alpha = 6 / R0
kMorses='1.6'

### setModel
#vdWs=(     'LJ' 'LJr8' 'LJr9' 'Morse' 'Buck') 
#vdW_flags=('1'  '2'    '3'    '4'     '5')
vdWs=(     'Morse')
vdW_flags=('4')

#elecs=(     'nil' 'Coul') 
#elec_flags=('0'   '1')
# 0=no Coul, 1=point charges, 2=soft clamping, 10-14=Boys clamping, 10=exact erf/r, 11=cubic C1, 12=quintic C2, 13=quartic even C1, 14=sextic even C2
icoul=1

#Hcorrs=(     'nil' 'H1' 'H2' 'H1H2')
#Hcorr_flags=('0'   '1'  '2'  '3')
#Hcorrs=(     'nil' 'H1' 'H2')
#Hcorr_flags=('0'   '1'  '2')
#Hcorrs=(     'nil')
#Hcorr_flags=('0')
#Hcorrs=(     'H1' 'H2' 'H1H2')
#Hcorr_flags=('1'  '2'  '3')
#Hcorrs=(     'H1' 'H2')
#Hcorr_flags=('1'  '2')
Hcorrs=(     'H1')
Hcorr_flags=('1')
#Hcorrs=(     'H2')
#Hcorr_flags=('2')

#Eps=(     'nil' 'SR' 'SR2')
#Ep_flags=('0'   '1'  '2')
#Eps=(     'nil')
#Ep_flags=('0')
#Eps=(     'SR' 'SR2')
#Ep_flags=('1'  '2')
Eps=(     'SR')
Ep_flags=('1')
#Eps=(     'SR2')
#Ep_flags=('2')

### LEpairs
#Lepairss='0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 2.0'  # distance of the Epair from the host atom
Lepairss='0.6'

### dof scans
dofnames=( 'H_O.H'             'E_H_O.R'           'E_H_O.H'     'O_3.H'              'E_O_3.R'     'E_O_3.H')
dofs=(     '0'                 '1'                 '2'           '3'                  '4'           '5')                  # degree of freedom to scan (see dofSelection file)
dofranges=('0.0 1.0 100'       '0.0 1.0 100'       '0.0 1.0 100' '-1.0 0.0 100'       '0.0 1.0 100' '-1.0 0.0 100')       # range of values of the degree of freedom to scan
dofstarts=('0.243309035386276' '0.246399924919899' '0.000001'    '-0.241077929836578' '0.000001'    '-0.223218554551512') # starting values of the degree of freedom to scan
ndofs=${#dofs[@]}

# setPenalty
regularize=0
regcountweight=0
softclamp=1
#softclamp_start=0.0 # eV ~ 0 kcal
#softclamp_max=0.0867 # eV ~ 2 kcal
#softclamp_start=0.0867 # eV ~ 2 kcal
#softclamp_max=0.1735 # eV ~ 4 kcal
softclamp_start=0.1735 # eV ~ 4 kcal
softclamp_max=0.2602 # eV ~ 6 kcal
user_weights=1
n_before=100
weight_a=1.0
weight_alpha=4.0
emin_min=-1.0 # eV
emin0=0.0 # eV

### folders
dir_main=/home/niko/work/HBOND
dir_exec=$dir_main/FireCore/tests/tFitREQ_PN
dir_run=$dir_exec/scan
dir_firecore=$dir_main/FireCore
dir_inputs=$dir_main/REFERENCE/2-pairs_small_small/4-to_firecore/confs_wb97m
dir_data=$dir_exec/data

### flags
# compile
rm_flag=false
#rm_flag=true
compile_flag=true
asan_flag=false
# run
run_flag=true
#run_flag=false
# plot
plot_flag=true
#show_flag=true
show_flag=false

#####################################################################################################################

### compile Firecore
if ( $compile_flag ) ; then
    cd $dir_firecore/cpp/Build/libs/Molecular
    if ( $rm_flag ) ; then rm -f libFitREQ_lib.so ; fi
    make FitREQ_PN_lib
    cd $dir_exec
    if ( $asan_flag ) ; then
	LD_PRELOAD=$(g++ -print-file-name=libasan.so)
	echo   $LD_PRELOAD
	export LD_PRELOAD
    fi
fi

for iv in ${!vdWs[@]} ; do
    ivdw=${vdW_flags[$iv]}
    vdw=${vdWs[$iv]}
    kMorses_tmp=$kMorses
    if [ ${vdw:0:2} == 'LJ' ] ; then kMorses_tmp=1.6 ; fi
    for kmorse in $kMorses_tmp ; do
	for ih in ${!Hcorrs[@]} ; do
	    ihbond=${Hcorr_flags[$ih]}
	    hbond=${Hcorrs[$ih]}
            for ie in ${!Eps[@]} ; do
		iepairs=${Ep_flags[$ie]}
		epairs=$iepairs
		if [ $iepairs -gt 0 ] ; then epairs=1 ; fi
		epairs_flag=${Eps[$ie]}
		#if [ $ihbond -eq 0 ] && [ $iepairs -eq 0 ] ; then continue ; fi
		Lepairss_tmp=$Lepairss
		if [ $iepairs -eq 0 ] ; then Lepairss_tmp=1.0 ; fi
		for lepairs in $Lepairss_tmp ; do
		    if [ $iepairs -eq 0 ] ; then
			if [ ${vdw:0:2} == 'LJ' ] ; then
			    longname="vdW$vdw-Hcorr$hbond-Ep$epairs_flag"
			else
			    longname="vdW$vdw-kMorse$kmorse-Hcorr$hbond-Ep$epairs_flag"
			fi
		    else
			if [ ${vdw:0:2} == 'LJ' ] ; then
			    longname="vdW$vdw-Hcorr$hbond-Ep$epairs_flag-LEpairs$lepairs"
			else
			    longname="vdW$vdw-kMorse$kmorse-Hcorr$hbond-Ep$epairs_flag-LEpairs$lepairs"
			fi
		    fi
		    echo ----------------------------------------------------
		    echo $longname
		    echo ----------------------------------------------------
		    
		    ### run
		    if ( $run_flag ) ; then
			# write dofSelection.dat files
			mkdir -p $dir_run
			rm -f $dir_run/dofSelection.dat
			touch $dir_run/dofSelection.dat
			for (( idof=0 ; idof < ndofs ; idof++ )) ; do
			    grep val_${dofnames[$idof]} $dir_data/dofSelection_run.dat | sed "s?val_${dofames[$idof]}?${dofstarts[$idof]}?" >> $dir_run/dofSelection.dat
			done
			# run
			mkdir -p $dir_run/dat
			for (( idof=0 ; idof < ndofs ; idof++ )) ; do
			    dof=${dofs[$idof]}	
			    python3 -u scan.py \
				    --ivdw $ivdw --icoul $icoul --ihbond $ihbond --epairs $epairs --iepairs $iepairs --lepairs $lepairs --kmorse $kmorse \
				    --dof_selection $dir_run/dofSelection.dat --inputs_dir $dir_inputs --inputs $files_inputs \
				    --regularize $regularize --regcountweight $regcountweight --softclamp $softclamp --softclamp_start $softclamp_start --softclamp_max $softclamp_max \
				    --user_weights $user_weights --n_before $n_before --weight_a $weight_a --weight_alpha $weight_alpha --emin_min $emin_min --emin0 $emin0 \
				    --mode scan --scan_dofs $dof --scan_range ${dofranges[$idof]} --scan_outfile $dir_run/DOF_scan 1> $dir_run/dat/run-$longname.out
			    mv $dir_run/DOF_scan_DOF${dof}_ana.dat $dir_run/dat/scan-$longname-${dofnames[$idof]}_ana.dat
			    mv $dir_run/DOF_scan_DOF${dof}_num.dat $dir_run/dat/scan-$longname-${dofnames[$idof]}_num.dat
#exit			    
			done
			rm $dir_run/dofSelection.dat
		    fi
		    
		    ### plots
		    if ( $plot_flag ) ; then
			mkdir -p $dir_run/png
			for (( idof=0 ; idof < ndofs ; idof++ )) ; do
			    dof=${dofs[$idof]}	
			    # fitness
			    echo "set terminal pngcairo noenhanced crop size 1400,900 font 'Times New Roman,18'" > $dir_run/x.gp
			    echo "set output '$dir_run/png/fitness-$longname-${dofnames[$idof]}.png'" >> $dir_run/x.gp
			    echo "set xlabel '${dofnames[$idof]} (a.u.)'" >> $dir_run/x.gp
			    echo "set ylabel 'Fitness (eV)'" >> $dir_run/x.gp
			    echo "plot '$dir_run/dat/scan-$longname-${dofnames[$idof]}_ana.dat' u 1:2 w l lw 3 lc rgb 'black' notitle" >> $dir_run/x.gp
			    gnuplot $dir_run/x.gp
			    if ( $show_flag ) ; then eog $dir_run/png/fitness-$longname-${dofnames[$idof]}.png & fi
			    # derivatives
			    echo "set terminal pngcairo noenhanced crop size 1400,900 font 'Times New Roman,18'" > $dir_run/x.gp
			    echo "set output '$dir_run/png/gradient-$longname-${dofnames[$idof]}.png'" >> $dir_run/x.gp
			    echo "set xlabel '${dofnames[$idof]} (a.u.)'" >> $dir_run/x.gp
			    echo "set ylabel 'Gradient (eV/a.u.)'" >> $dir_run/x.gp
			    echo "plot '$dir_run/dat/scan-$longname-${dofnames[$idof]}_ana.dat' u 1:3 w p pt 2 ps 3 lc rgb 'black' title 'Analytical',"'\' >> $dir_run/x.gp
			    echo "     '$dir_run/dat/scan-$longname-${dofnames[$idof]}_num.dat' u 1:2 w l lw 3      lc rgb 'black' title 'Numerical'" >> $dir_run/x.gp
			    gnuplot $dir_run/x.gp
			    if ( $show_flag ) ; then eog $dir_run/png/gradient-$longname-${dofnames[$idof]}.png & fi
			    rm $dir_run/x.gp
			done
		    fi
		    
		done
	    done
	done
    done
done

exit
