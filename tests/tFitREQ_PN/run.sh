#!/bin/bash

### inputs
files_inputs='H2O-A1_H2O-D1-y.xyz'
#files_inputs='H2O-A1_H2O-D1-y.xyz H2O-A1_H2O-D1-z.xyz H2O-D1_H2O-A1-y.xyz H2O-D1_H2O-A1-z.xyz'
#files_inputs='CH2O-A1_NH3-D1-y.xyz CH2O-A1_NH3-D1-z.xyz NH3-D1_CH2O-A1-y.xyz NH3-D1_CH2O-A1-z.xyz NH3-A1_NH3-D1-y.xyz NH3-A1_NH3-D1-z.xyz NH3-D1_NH3-A1-y.xyz NH3-D1_NH3-A1-z.xyz'

### kMorse
#kMorses='1.6 1.7 1.8 -1' # value of alpha for Morse and Buck interactions, a negative value means that alpha = 6 / R0
#kMorses='1.7 1.8 -1' 
kMorses='1.7' 

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

#Eps=(     'nil' 'SR' 'SR2' 'SR3')
#Ep_flags=('0'   '1'  '2'   '3')
#Eps=(     'nil')
#Ep_flags=('0')
Eps=(     'SR' 'SR2' 'SR3')
Ep_flags=('1'  '2'   '3')
#Eps=(     'SR')
#Ep_flags=('1')
#Eps=(     'SR2')
#Ep_flags=('2')
#Eps=(     'SR3')
#Ep_flags=('3')

### LEpairs
#Lepairss='0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 2.0'  # distance of the Epair from the host atom
#Lepairss='0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 2.0'
Lepairss='0.5 1.0 1.1 1.2'
#Lepairss='0.6'

### initial values for fitting
dofnames=( 'H_O.H' 'E_H_O.R' 'E_H_O.H' 'O_3.H' 'E_O_3.R' 'E_O_3.H')
dofstarts=('0.5'   '0.5'     '1.0'     '-0.5'  '0.5'     '-1.0')
##dofstarts=('0.52'   '0.99'     '1.71'     '-1.31'  '1.31'     '-3.62')
#dofnames=( 'H_N.H' 'E_H_N.R' 'E_H_N.H' 'N_3.H' 'E_N_3.R' 'E_N_3.H' 'O_2.H' 'E_O_2.R' 'E_O_2.H')
#dofstarts=('0.5'   '0.5'     '1.0'     '-0.5'  '0.5'     '-1.0'    '-0.5'  '0.5'     '-1.0')

### energy for plotting
DEmax=1.0 # kcal/mol, for plotting: Emax = Emin + DEmax if DEmax > 0.0, Emax = -Emin if DEmax < 0.0
#DEmax=-1

### setPenalty
clamp=1 # Hardly restrain the values of parameters during optimization, 0=no, 1=yes
regularize=0
regcountweight=0
softclamp=1
#softclamp_start=0.0 # eV ~ 0 kcal
#softclamp_max=0.0867 # eV ~ 2 kcal
softclamp_start=0.0867 # eV ~ 2 kcal
softclamp_max=0.1735 # eV ~ 4 kcal
#softclamp_start=0.1735 # eV ~ 4 kcal
#softclamp_max=0.2602 # eV ~ 6 kcal
user_weights=1
n_before=100
weight_a=1.0
weight_alpha=4.0
emin_min=-1.0 # eV
emin0=0.0 # eV

### atom types
atom_types_nofit='H_a C_R C_2 C_1'
atom_types_all=('H_C1' 'H_N' 'H_O' 'H_F' 'H_Br' 'H_Cl' 'N_3' 'N_R' 'N_2' 'N_1' 'O_3' 'O_2' 'F_'  'Br'  'Cl')
fits=(          false  false true  false false  false  false false false false true  false false false false) # H2O
#fits=(          false  true  false false false  false  true  false false false false true  false false false) # NH3 CH2O
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
#plot_flag=false
#show_flag=true
show_flag=false

### optimization algorithm
#ialgs='0 1 2 3' # 0=move_GD, 1=move_MD, 2=move_GD_BB_short, 3=move_GD_BB_long
ialg='2'

### distances for plot
dmin=1.4 # Ang
dmax=3.0 # Ang

### folders
dir_main=/home/niko/work/HBOND
dir_exec=$dir_main/FireCore/tests/tFitREQ_PN
dir_run=$dir_exec/run
#dir_mols=$dir_main/REFERENCE/1-small_molecules/1-from_avogadro
dir_firecore=$dir_main/FireCore
dir_inputs=$dir_main/REFERENCE/2-pairs_small_small/4-to_firecore/confs_wb97m
dir_data=$dir_exec/data

### molecules
#mols='C4H3NO2 C4H5N C5H5N CH2NH CH2O H2O HCN HCONH2 HCOOH HF NH3' # HBr HCl

### python options
# setOptimization
nstep=2000
#nstep=100000
#fmax=1.0e-8
#dt=0.02
#max_step=0.05
#damping=0.01
#iparallel=1 # 0=serial, 1=serial function with bOMP=1
# setOutput
#outxyz=0
#savejustelementxyz=0
#outxyz_fname=$dir_run/input_epairs.xyz

#####################################################################################################################

n_files_inputs=0 ; for f in $files_inputs ; do n_files_inputs=$(( n_files_inputs + 1 )) ; done
    
mkdir -p $dir_run
gfortran -Wall $dir_exec/minmax.f90 -o $dir_run/minmax.x
gfortran -Wall $dir_exec/minmaxdiff.f90 -o $dir_run/minmaxdiff.x
gfortran -Wall $dir_exec/rmse.f90 -o $dir_run/rmse.x
source utils.sh

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

# count runs
count_runs

irun=0
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
		    irun=$(( irun + 1 ))
		    echo ----------------------------------------------------
		    echo $irun/$nruntot $longname
		    echo ----------------------------------------------------
			
		    ### atom types
		    atom_types_begin
		    
		    ### autoinputs # TO BE CHECKED
		    # ALSO THERE CAN BE TWO MODES TO RUN:
		    # 1) SELECT SOME ATOM TYPES AND AUTOMATICALLY INCLUDE ALL RELEVANT MOLECULES
		    # 2) SELECT SOME MOLECULE AND AUTOMATICALLY INCLUDE ALL INVOLVED ATOM TYPES
		    #autoinputs

		    ### run
		    if ( $run_flag ) ; then
			check_tbd_run $dir_run/dat/run-$longname.out
			if ( $tbd ) ; then
			    # write dofSelection.dat files
			    rm -f $dir_run/dofSelection.dat
			    touch $dir_run/dofSelection.dat
			    for ipar in ${!startpars[@]} ; do
				grep " val_${parnames[$ipar]} " $dir_data/dofSelection_run.dat | sed "s?val_${parnames[$ipar]}?${startpars[$ipar]}?" >> $dir_run/dofSelection.dat
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
			    print_output
			    # output store
			    store_output
			fi
		    fi
		    
		    ### plots
		    if ( $plot_flag ) ; then
			mkdir -p $dir_run/png
			# traj
			if [ ${#startpars[@]} -gt 0 ] ; then
			    check_tbd_plot $dir_run/png/opt-$longname.png
			    if ( $tbd ) ; then
				echo "  plot traj $longname"
				plot_traj 
				if ( $show_flag ) ; then eog $dir_run/png/opt-$longname.png & fi
			    fi
			fi
			# polar energy maps
			check_tbd_plot $dir_run/png/polarmap-$longname.png
			if ( $tbd ) ; then
			    echo "  plot polarmap $longname"
			    for f1 in $files_inputs ; do
				f=${f1//.xyz}
				convert -background White -fill black -font Helvetica -pointsize 36 -gravity center label:"$f" $dir_run/title.png
				find_eminmax $dir_run/dat/${f}__ref-$longname.dat $DEmax
				plot_polarmap ref ref ${f}__ref
				plot_polarmap model model ${f}__model
				if (( $(echo "$DEmax > 0.0" | bc -l) )); then
				    Emin=-$DEmax
				    Emax=$DEmax
				else
				    find_eminmaxdiff $dir_run/dat/${f}__diff-$longname.dat
				fi
				calc_rmse $dir_run/dat/${f}__diff-$longname.dat
				plot_polarmap diff "RMSE = $rmse kcal/mol" ${f}__diff
				convert -background white -gravity North $dir_run/polar_ref.png $dir_run/polar_model.png $dir_run/polar_diff.png +append $dir_run/png/polarmap_$f-$longname.png
				convert $dir_run/png/polarmap_$f-$longname.png -background White -splice 0x60 $dir_run/input_spaced.png
				convert $dir_run/input_spaced.png $dir_run/title.png -gravity north -geometry +0+0 -composite $dir_run/polarmap_$f.png
				rm $dir_run/title.png $dir_run/polar_ref.png $dir_run/polar_model.png $dir_run/polar_diff.png $dir_run/png/polarmap_$f-$longname.png $dir_run/input_spaced.png
			    done
			    if [ $n_files_inputs -eq 1 ] ; then
				mv $dir_run/polarmap_$f.png $dir_run/png/polarmap-$longname.png
			    else
				com='convert -background white -gravity North'
				for f1 in $files_inputs ; do
				    com="$com $dir_run/polarmap_${f1//.xyz}.png"
				done
				com="$com -append $dir_run/png/polarmap-$longname.png"
				eval "$com"
				rm $dir_run/polarmap_*.png
			    fi
			    if ( $show_flag ) ; then eog $dir_run/png/polarmap-$longname.png & fi
			fi
			# cartesian energy maps
			check_tbd_plot $dir_run/png/cartmap-$longname.png
			if ( $tbd ) ; then
			    echo "  plot cartmap $longname"
			    for f1 in $files_inputs ; do
				f=${f1//.xyz}
				convert -background White -fill black -font Helvetica -pointsize 36 -gravity center label:"$f" $dir_run/title.png
				find_eminmax $dir_run/dat/${f}__ref-$longname.dat $DEmax
				plot_cartmap  ref ref ${f}__ref
				plot_cartmap  model model ${f}__model
				if (( $(echo "$DEmax > 0.0" | bc -l) )); then
				    Emin=-$DEmax
				    Emax=$DEmax
				else
				    find_eminmaxdiff $dir_run/dat/${f}__diff-$longname.dat
				fi
				calc_rmse $dir_run/dat/${f}__diff-$longname.dat
				plot_cartmap  diff "RMSE = $rmse kcal/mol" ${f}__diff
				convert -background white -gravity North $dir_run/cart_ref.png $dir_run/cart_model.png $dir_run/cart_diff.png +append $dir_run/png/cartmap_$f-$longname.png
				convert $dir_run/png/cartmap_$f-$longname.png -background White -splice 0x60 $dir_run/input_spaced.png
				convert $dir_run/input_spaced.png $dir_run/title.png -gravity north -geometry +0+0 -composite $dir_run/cartmap_$f.png
				rm $dir_run/title.png $dir_run/cart_ref.png $dir_run/cart_model.png $dir_run/cart_diff.png $dir_run/png/cartmap_$f-$longname.png $dir_run/input_spaced.png
			    done
			    if [ $n_files_inputs -eq 1 ] ; then
				mv $dir_run/cartmap_$f.png $dir_run/png/cartmap-$longname.png
			    else
				com='convert -background white -gravity North'
				for f1 in $files_inputs ; do
				    com="$com $dir_run/cartmap_${f1//.xyz}.png"
				done
				com="$com -append $dir_run/png/cartmap-$longname.png"
				eval "$com"
				rm $dir_run/cartmap_*.png
			    fi
			    if ( $show_flag ) ; then eog $dir_run/png/cartmap-$longname.png & fi
			fi
			# min lines
			check_tbd_plot $dir_run/png/min-$longname.png 
			if ( $tbd ) ; then
			    echo "  plot minlines $longname"
			    for f1 in $files_inputs ; do
				f=${f1//.xyz}
				convert -background White -fill black -font Helvetica -pointsize 36 -gravity center label:"$f" $dir_run/title.png
				plot_minlines
				convert $dir_run/png/min_$f-$longname.png -background White -splice 0x60 $dir_run/input_spaced.png
				convert $dir_run/input_spaced.png $dir_run/title.png -gravity north -geometry +0+0 -composite $dir_run/min_$f.png
				rm $dir_run/title.png $dir_run/png/min_$f-$longname.png $dir_run/input_spaced.png
			    done
			    if [ $n_files_inputs -eq 1 ] ; then
				mv $dir_run/min_$f.png $dir_run/png/min-$longname.png
			    else
				com='convert -background white -gravity North'
				for f1 in $files_inputs ; do
				    com="$com $dir_run/min_${f1//.xyz}.png"
				done
				com="$com -append $dir_run/png/min-$longname.png"
				eval "$com"
				rm $dir_run/min_*.png
			    fi
			    if ( $show_flag ) ; then eog $dir_run/png/min-$longname.png & fi
			fi
		    fi
#exit
		    atom_types_end
		done
	    done
	done
    done
done

rm $dir_run/minmax.x $dir_run/minmaxdiff.x $dir_run/rmse.x

exit
