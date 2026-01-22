#!/bin/bash

### inputs
files_inputs='H2O-A1_H2O-D1-y.xyz'
#files_inputs='H2O-A1_H2O-D1-y.xyz H2O-A1_H2O-D1-z.xyz H2O-D1_H2O-A1-y.xyz H2O-D1_H2O-A1-z.xyz'

### kMorse
kMorses='1.6 -1' # value of alpha for Morse and Buck interactions, a negative value means that alpha = 6 / R0
#kMorses='1.6' 
#kMorses='-1' 

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
Hcorrs=(     'H1' 'H2')
Hcorr_flags=('1'  '2')
#Hcorrs=(     'H1')
#Hcorr_flags=('1')
#Hcorrs=(     'H2')
#Hcorr_flags=('2')

#Eps=(     'nil' 'SR' 'SR2')
#Ep_flags=('0'   '1'  '2')
#Eps=(     'nil')
#Ep_flags=('0')
Eps=(     'SR' 'SR2')
Ep_flags=('1'  '2')
#Eps=(     'SR')
#Ep_flags=('1')
#Eps=(     'SR2')
#Ep_flags=('2')

### LEpairs
Lepairss='0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 2.0'  # distance of the Epair from the host atom
#Lepairss='1.0 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 2.0'
#Lepairss='0.5 1.0 1.5'
#Lepairss='1.6'
#Lepairss='1.6 1.7 1.8 1.9 2.0'

### initial values for fitting (the same for all atoms)
Hat=0.5
Re=0.5
He=1.0

# setPenalty
clamp=1 # Hardly restrain the values of parameters during optimization, 0=no, 1=yes
regularize=0
regcountweight=0
softclamp=1
#softclamp_start=0.0 # 0 kcal
#softclamp_max=0.0867 # 2 kcal
#softclamp_start=0.0867 # 2 kcal
#softclamp_max=0.1735 # 4 kcal
softclamp_start=0.1735 # 4 kcal
softclamp_max=0.2602 # 6 kcal
user_weights=1
n_before=100
weight_a=1.0
weight_alpha=4.0
emin_min=-1
emin0=0.0

### optimization algorithm
#ialgs='0 1 2 3' # 0=move_GD, 1=move_MD, 2=move_GD_BB_short, 3=move_GD_BB_long
ialg='2'

### distances for plot
dmin=1.4 # Ang
dmax=3.0 # Ang

### folders
dir_main=/home/niko/work/HBOND
dir_run=$dir_main/FIT/run
#dir_mols=$dir_main/REFERENCE/1-small_molecules/1-from_avogadro
dir_firecore=$dir_main/FireCore
dir_inputs=$dir_main/REFERENCE/2-pairs_small_small/4-to_firecore/confs_wb97m
dir_data=$dir_main/FIT/data

### molecules
#mols='C4H3NO2 C4H5N C5H5N CH2NH CH2O H2O HCN HCONH2 HCOOH HF NH3' # HBr HCl

### atom types
atom_types_nofit='H_a C_R C_2 C_1'
atom_types_all=('H_C1' 'H_N' 'H_O' 'H_F' 'H_Br' 'H_Cl' 'N_3' 'N_R' 'N_2' 'N_1' 'O_3' 'O_2' 'F_'  'Br'  'Cl')
fits=(          false  false true  false false  false  false false false false true  false false false false) # H2O
natom_types_all=${#atom_types_all[@]}

### flags
# compile
rm_flag=false
#rm_flag=true
compile_flag=true
asan_flag=false
# autoinputs
#autoinputs_flag=true
#autoinputs_flag=false
# run
run_flag=true
#run_flag=false
# plot
plot_flag=true
#show_flag=true
show_flag=false

### python options
# setOptimization
nstep=2000
#fmax=1.0e-8
#dt=0.02
#max_step=0.05
#damping=0.01
#iparallel=0 # 0=serial, 1=serial function with bOMP=1
# setOutput
#outxyz=0
#savejustelementxyz=0
#outxyz_fname=$dir_run/input_epairs.xyz

#####################################################################################################################

gfortran -Wall minmax.f90 -o minmax.x
source utils.sh

### compile Firecore
if ( $compile_flag ) ; then
    cd $dir_firecore/cpp/Build/libs/Molecular
    if ( $rm_flag ) ; then rm -f libFitREQ_lib.so ; fi
    make FitREQ_lib
    cd - 1> /dev/null
    if ( $asan_flag ) ; then
	LD_PRELOAD=$(g++ -print-file-name=libasan.so)
	echo   $LD_PRELOAD
	export LD_PRELOAD
    fi
fi

for kmorse in $kMorses ; do
    for iv in ${!vdWs[@]} ; do
	ivdw=${vdW_flags[$iv]}
	vdw=${vdWs[$iv]}
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
		    longname="kMorse$kmorse-vdW$vdw-Hcorr$hbond-Ep$epairs_flag-LEpairs$lepairs"
		    echo ------------------------------------------
		    echo $longname
		    echo ------------------------------------------
			
		    ### atom types
		    declare -a atom_types
		    declare -a parnames
		    declare -a startpars
		    declare -a parunits
		    # consider just atom types to be fitted
		    atom_types_list=''
		    for (( ia=0 ; ia < natom_types_all ; ia++ )) ; do
			at=${atom_types_all[$ia]} ;
			fit=${fits[$ia]} ;
			if ( ! $fit ) ; then continue ; fi
			atom_types=( "${atom_types[@]}" "$at" )
			atom_types_list="$atom_types_list $at"
			if [ $ihbond -gt 0 ] ; then
			    parnames=( "${parnames[@]}" "${at}.H" )
			    if [ "${at:0:1}" == 'H' ] ; then
				startpars=( "${startpars[@]}" "$Hat" )
			    else
				startpars=( "${startpars[@]}" "-$Hat" )
			    fi
			    parunits=( "${parunits[@]}" "-" )
			fi
			if [ $iepairs -gt 0 ] ; then 
			    parnames=( "${parnames[@]}" "E_${at}.R" "E_${at}.H" )
			    if [ "${at:0:1}" == 'H' ] ; then
				startpars=( "${startpars[@]}" "$Re" "$He" )
			    else
				startpars=( "${startpars[@]}" "$Re" "-$He" )
			    fi
			    parunits=( "${parunits[@]}" "Ang" "eV^1/2" )
			fi
		    done
		    natom_types=${#atom_types[@]}
		    npars=${#startpars[@]}
		    
		    ### autoinputs # TO BE CHECKED
		    # ALSO THERE CAN BE TWO MODES TO RUN:
		    # 1) SELECT SOME ATOM TYPES AND AUTOMATICALLY INCLUDE ALL RELEVANT MOLECULES
		    # 2) SELECT SOME MOLECULE AND AUTOMATICALLY INCLUDE ALL INVOLVED ATOM TYPES
		    # find molecules
		    #mols_to_fit=''
		    #for mol in $mols ; do
		    #	types=`head -2 $dir_mols/$mol.xyz | tail -1`
		    #	for a in $atom_types_nofit ; do types=`echo $types | sed "s?$a??g"` ; done
		    #	for at in "${atom_types[@]}" ; do
		    #	    types=`echo $types | sed "s?$at??g"`
		    #	done
		    #	if [ "${types// }" == '' ] ; then mols_to_fit="$mols_to_fit $mol" ; fi
		    #done
		    ## find pairs
		    #files_inputs=''
		    #for pair_xyz in `ls $dir_inputs/*.xyz | sed "s?$dir_inputs/??g"` ; do
		    #	pair=`echo $pair_xyz | sed "s?-? ?g" | sed "s?_? ?g" | awk '{print $1,$3}'`
		    #	for mol in $mols_to_fit ; do
		    #	    pair=`echo $pair | sed "s?$mol??g"`
		    #	done
		    #	if [ "${pair// }" == '' ] ; then files_inputs="$files_inputs $pair_xyz" ; fi
		    #done
		    
		    ### run
		    if ( $run_flag ) ; then
			# write dofSelection.dat files
			mkdir -p $dir_run
			rm -f $dir_run/dofSelection.dat
			for ipar in ${!startpars[@]} ; do
			    grep val_${parnames[$ipar]} $dir_data/dofSelection_template.dat.sed | sed "s?val_${parnames[$ipar]}?${startpars[$ipar]}?" >> $dir_run/dofSelection.dat
			done
			# run
			mkdir -p $dir_run/dat
			python3 -u run.py \
				--ivdw $ivdw --icoul $icoul --ihbond $ihbond --epairs $epairs --iepairs $iepairs --lepairs $lepairs --kmorse $kmorse \
				--dof_selection $dir_run/dofSelection.dat --inputs_dir $dir_inputs --inputs $files_inputs \
				--clamp $clamp --regularize $regularize --regcountweight $regcountweight \
				--softclamp $softclamp --softclamp_start $softclamp_start --softclamp_max $softclamp_max \
				--user_weights $user_weights --n_before $n_before --weight_a $weight_a --weight_alpha $weight_alpha --emin_min $emin_min --emin0 $emin0 \
				--ialg $ialg --nstep $nstep --out_dir $dir_run 1> $dir_run/dat/run-$longname.out
#exit			    
			# output print
			print_output $dir_run $longname $ihbond $iepairs
			# output store
			for ipar in ${!startpars[@]} ; do
			    mv $dir_run/${parnames[$ipar]}.dat $dir_run/dat/${parnames[$ipar]}-$longname.dat
			done
			for f1 in $files_inputs ; do
			    f=${f1//.xyz}
			    mv $dir_run/${f}__diff.dat        $dir_run/dat/${f}__diff-$longname.dat
			    mv $dir_run/${f}__model.dat       $dir_run/dat/${f}__model-$longname.dat
			    mv $dir_run/${f}__model_lines.dat $dir_run/dat/${f}__model_lines-$longname.dat
			    mv $dir_run/${f}__ref.dat         $dir_run/dat/${f}__ref-$longname.dat
			    mv $dir_run/${f}__ref_lines.dat   $dir_run/dat/${f}__ref_lines-$longname.dat
			done
			rm $dir_run/dofSelection.dat
		    fi
		    
		    ### plots
		    if ( $plot_flag ) ; then
			mkdir -p $dir_run/png
			# paramater traj
			echo "set terminal pngcairo noenhanced crop size 1150,900 font 'Arial,18'" > x.gp
			echo "set output '$dir_run/pars.png'" >> x.gp
			echo "set title 'Optimization trajectory'" >> x.gp
			echo "set xlabel 'Steps'" >> x.gp
			echo "set ylabel 'Parameter value'" >> x.gp
			echo "plot "'\' >> x.gp
			for (( ipar=0 ; ipar < npars-1 ; ipar++ )) ; do
			    echo "'$dir_run/dat/${parnames[$ipar]}-$longname.dat'    u 1:2 w l lw 3 dt 1 lc $((ipar+1))  title '${parnames[$ipar]}',"'\' >> x.gp
			done
			echo "    '$dir_run/dat/${parnames[$npars-1]}-$longname.dat' u 1:2 w l lw 3 dt 1 lc $npars  title '${parnames[$npars-1]}'" >> x.gp
			gnuplot x.gp
			# traj
			echo "set terminal pngcairo noenhanced crop size 1400,900 font 'Arial,18'" > x.gp
			echo "set output '$dir_run/traj.png'" >> x.gp
			echo "set xlabel 'Steps'" >> x.gp
			echo "set ylabel 'Fitness (eV^2)' tc rgb 'black'" >> x.gp
			echo "set ytics nomirror tc rgb 'black'" >> x.gp
			echo "set logscale y" >> x.gp
			echo "set y2label 'Force norm (eV^2/a.u.^2)' tc rgb 'red'" >> x.gp
			echo "set y2tics nomirror tc rgb 'red'" >> x.gp
			echo "set logscale y2" >> x.gp
			echo "plot '$dir_run/dat/traj-$longname.dat' u 1:3 w l lw 3 lc rgb 'red'   axes x1y2 notitle,"'\' >> x.gp
			echo "     '$dir_run/dat/traj-$longname.dat' u 1:2 w l lw 3 lc rgb 'black' axes x1y1 notitle" >> x.gp
			gnuplot x.gp
			convert -background white -gravity North $dir_run/pars.png $dir_run/traj.png -append $dir_run/png/opt-$longname.png
			if ( $show_flag ) ; then eog $dir_run/png/opt-$longname.png & fi
			# energy maps
			for f1 in $files_inputs ; do
			    f=${f1//.xyz}
			    # ref polar
			    Emin=`grep -v nan $dir_run/dat/${f}__ref-$longname.dat | grep -v angle | awk 'BEGIN{a=1e99}{if ($3+0<a+0) a=$3} END{print a}'`
			    #Emax=$(awk -v emin="$Emin" 'BEGIN{print -emin}')
			    Emax=$(awk -v emin="$Emin" 'BEGIN{print emin+1.0}')
			    plot_polarmap $dir_run ref DFT $dmin $dmax $Emin $Emax ${f}__ref $longname
			    # ref cart
			    grep -v nan $dir_run/dat/${f}__ref-$longname.dat | grep -v angle > x.0
			    echo "x.0 x.1 $dmin $dmax" | ./minmax.x > /dev/null
			    sed "s?MYFILENAME?x.1?" sp_plot_polar.py | sed "s?MYMIN?$Emin?" | sed "s?MYMAX?$Emax?" | \
                                sed "s?MYTITLE?ref?" | sed "s?MYPNG?$dir_run/cart_ref.png?" > x.py
                            python x.py
			    # model polar
			    plot_polarmap $dir_run model $longname $dmin $dmax $Emin $Emax ${f}__model $longname
			    # model cart
			    grep -v nan $dir_run/dat/${f}__model-$longname.dat | grep -v angle > x.0
			    echo "x.0 x.1 $dmin $dmax" | ./minmax.x > /dev/null
			    sed "s?MYFILENAME?x.1?" sp_plot_polar.py | sed "s?MYMIN?$Emin?" | sed "s?MYMAX?$Emax?" | \
                                sed "s?MYTITLE?$longname?" | sed "s?MYPNG?$dir_run/cart_model.png?" > x.py
                            python x.py
			    # diff polar
			    Emin=`grep -v nan $dir_run/dat/${f}__diff-$longname.dat | grep -v angle | awk 'BEGIN{a=1e99}{if ($3+0<a+0) a=$3} END{print a}'`
			    Emax=`grep -v nan $dir_run/dat/${f}__diff-$longname.dat | grep -v angle | awk 'BEGIN{a=-1e99}{if ($3+0>a+0) a=$3} END{print a}'`
			    awk -v emin="$Emin" -v emax="$Emax" 'BEGIN{if (!(emin<0 && emax>0)) exit 1}' || { echo "ERROR: Emin must be <0 and Emax >0 (got Emin=$Emin, Emax=$Emax)"; exit 1; }
			    if awk -v emin="$Emin" -v emax="$Emax" 'BEGIN{exit !(-emin > emax)}'; then
				Emax=$(awk -v emin="$Emin" 'BEGIN{print -emin}')
			    else
				Emin=$(awk -v emax="$Emax" 'BEGIN{print -emax}')
			    fi
			    rmse=`grep -v nan $dir_run/dat/${f}__diff-$longname.dat | grep -v angle | awk '{a+=$3^2 ; printf "%25.2f\n", (a/NR)^0.5}' | tail -1 | awk '{print $1}'`
			    plot_polarmap $dir_run diff "RMSE = $rmse kcal/mol" $dmin $dmax $Emin $Emax ${f}__diff $longname
			    # diff cart
			    grep -v nan $dir_run/dat/${f}__diff-$longname.dat | grep -v angle > x.0
			    echo "x.0 x.1 $dmin $dmax" | ./minmax.x > /dev/null
			    sed "s?MYFILENAME?x.1?" sp_plot_polar.py | sed "s?MYMIN?$Emin?" | sed "s?MYMAX?$Emax?" | \
                                sed "s?MYTITLE?RMSE = $rmse kcal/mol?" | sed "s?MYPNG?$dir_run/cart_diff.png?" > x.py
                            python x.py
			    convert -background white -gravity North $dir_run/ref.png $dir_run/model.png $dir_run/diff.png +append $dir_run/png/polarmap_$f-$longname.png
			    if ( $show_flag ) ; then eog $dir_run/png/map_$f-$longname.png & fi
			    convert -background white -gravity North $dir_run/cart_ref.png $dir_run/cart_model.png $dir_run/cart_diff.png +append $dir_run/png/cartmap_$f-$longname.png
			    if ( $show_flag ) ; then eog $dir_run/png/cartmap_$f-$longname.png & fi
			done
			# min lines
			for f1 in $files_inputs ; do
			    f=${f1//.xyz}
			    # rmin
			    echo "set terminal pngcairo noenhanced crop size 700,900 font 'Arial,18'" > x.gp
			    echo "set output '$dir_run/rmin.png'" >> x.gp
			    echo "set xlabel 'Angle (deg)'" >> x.gp
			    echo "set xrange [-90:90]" >> x.gp
			    echo "set xtics 30" >> x.gp
			    echo "set ylabel 'Distance of the minimum (Ang)'" >> x.gp
			    echo "set yrange [$dmin:$dmax]" >> x.gp
			    echo "set logscale y" >> x.gp
			    echo "set ytics 1 nologscale" >> x.gp
			    echo "plot '$dir_run/dat/${f}__ref_lines-$longname.dat'   u 1:2 w lp lw 3 pt 6 ps 2 lc rgb 'black' title 'ref',"'\' >> x.gp
			    echo "     '$dir_run/dat/${f}__model_lines-$longname.dat' u 1:2 w lp lw 3 pt 6 ps 2 lc rgb 'red'   title 'model'" >> x.gp
			    gnuplot x.gp
			    # Emin
			    echo "set terminal pngcairo noenhanced crop size 700,900 font 'Arial,18'" > x.gp
			    echo "set output '$dir_run/emin.png'" >> x.gp
			    echo "set xlabel 'Angle (deg)'" >> x.gp
			    echo "set xrange [-90:90]" >> x.gp
			    echo "set xtics 30" >> x.gp
			    echo "set ylabel 'Energy of the minimum (kcal/mol)'" >> x.gp
			    echo "plot '$dir_run/dat/${f}__ref_lines-$longname.dat'   u 1:3 w lp lw 3 pt 6 ps 2 lc rgb 'black' title 'ref',"'\' >> x.gp
			    echo "     '$dir_run/dat/${f}__model_lines-$longname.dat' u 1:3 w lp lw 3 pt 6 ps 2 lc rgb 'red'   title 'model'" >> x.gp
			    gnuplot x.gp
			    convert -background white -gravity North $dir_run/rmin.png $dir_run/emin.png +append $dir_run/png/min_$f-$longname.png
			    if ( $show_flag ) ; then eog $dir_run/png/min_$f-$longname.png & fi
			done
			rm x.gp $dir_run/pars.png $dir_run/traj.png $dir_run/ref.png $dir_run/model.png $dir_run/diff.png $dir_run/rmin.png $dir_run/emin.png
			rm x.0 x.1 x.py $dir_run/cart_ref.png $dir_run/cart_model.png $dir_run/cart_diff.png 
		    fi
		    
		    unset atom_types
		    unset parnames
		    unset startpars
		    unset parunits
		done
	    done
	done
    done
done

exit
