#!/bin/bash

#----------------------------------------------------
#444/1419 vdWMorse-kMorse1.6-HcorrH1-EpSR-LEpairs1.2
#----------------------------------------------------
#CONVERGED in 273 iterations (|F|=9.58922e-09 < F2max= 1e-08) 
#VERY FINAL |E|=0.000114507004453663  DOFs= 0.300902933264142 0.019520964588372 0.724021869024964 -0.571304074631488 1.62981149989273 -1.22756624272456 
#VERY FINAL |F|=9.58921546251403e-09 fDOFs= 1.64089440768538e-11 -9.57650777722579e-09 -2.1091786003435e-11 -2.43064367179314e-11 -4.27050564018385e-10 2.44693431205468e-10 
#DOF   0  t:  38  H_O      H :               0.300902933264142     1.64e-11 
#DOF   1  t:  12  E_H_O    R :               0.019520964588372    -9.58e-09 
#DOF   2  t:  12  E_H_O    H :               0.724021869024964    -2.11e-11 
#DOF   3  t:  30  O_3      H :              -0.571304074631488    -2.43e-11 
#DOF   4  t:   6  E_O_3    R :               1.629811499892730    -4.27e-10 
#DOF   5  t:   6  E_O_3    H :              -1.227566242724565     2.45e-10 
#----------------------------------------------------
#445/1419 vdWMorse-kMorse1.6-HcorrH1-EpSR-LEpairs1.3
#----------------------------------------------------
#step= 1999 dt= 0.821172
#VERY FINAL |E|=0.000176636505470402  DOFs= 1e-10 1.13505831999308 0.895637289236092 -0.166770952165735 0.182566560262891 -0.907237476833397 
#VERY FINAL |F|=0.00114004448602431 fDOFs= -0.00114004225207844 8.91575600423956e-08 8.5062220058451e-08 -2.20561701268363e-07 -2.24229168959049e-06 4.34304745660746e-08 
#DOF   0  t:  38  H_O      H :               0.000000000100000    -1.14e-03 
#DOF   1  t:  12  E_H_O    R :               1.135058319993082     8.92e-08 
#DOF   2  t:  12  E_H_O    H :               0.895637289236092     8.51e-08 
#DOF   3  t:  30  O_3      H :              -0.166770952165735    -2.21e-07 
#DOF   4  t:   6  E_O_3    R :               0.182566560262891    -2.24e-06 
#DOF   5  t:   6  E_O_3    H :              -0.907237476833397     4.34e-08 

### LEpairs
Lepairss=('1.2' '1.3')
nLepairss=${#Lepairss[@]}

### actual values of parameters
dofnames=('H_O.H'             'E_H_O.R'           'E_H_O.H'           'O_3.H'              'E_O_3.R'           'E_O_3.H')
dofs1=(   '0.300902933264142' '0.019520964588372' '0.724021869024964' '-0.571304074631488' '1.629811499892730' '-1.227566242724565')
dofs2=(   '0.000000000100000' '1.135058319993082' '0.895637289236092' '-0.166770952165735' '0.182566560262891' '-0.907237476833397')
ndofs=${#dofs[@]}
nsteps=100

### inputs
files_inputs='H2O-A1_H2O-D1-y.xyz'
#files_inputs='H2O-A1_H2O-D1-y.xyz H2O-A1_H2O-D1-z.xyz H2O-D1_H2O-A1-y.xyz H2O-D1_H2O-A1-z.xyz'

### kMorse
#kMorses='1.6 1.7 1.8 -1' # value of alpha for Morse and Buck interactions, a negative value means that alpha = 6 / R0
kmorse='1.6' 

### setModel
#vdWs=(     'LJ' 'LJr8' 'LJr9' 'Morse' 'Buck') 
#vdW_flags=('1'  '2'    '3'    '4'     '5')
ivdw=4
vdw='Morse'

#elecs=(     'nil' 'Coul') 
#elec_flags=('0'   '1')
# 0=no Coul, 1=point charges, 2=soft clamping, 10-14=Boys clamping, 10=exact erf/r, 11=cubic C1, 12=quintic C2, 13=quartic even C1, 14=sextic even C2
icoul=1

#Hcorrs=(     'nil' 'H1' 'H2' 'H1H2')
#Hcorr_flags=('0'   '1'  '2'  '3')
ihbond=1
hbond='H1'

#Eps=(     'nil' 'SR' 'SR2')
#Ep_flags=('0'   '1'  '2')
iepairs=1
epairs_flag='SR'
epairs=1

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
dir_run=$dir_exec/rigidscan
dir_firecore=$dir_main/FireCore
dir_inputs=$dir_main/REFERENCE/2-pairs_small_small/4-to_firecore/confs_wb97m
dir_data=$dir_exec/data

### flags
# compile
rm_flag=false
#rm_flag=true
#compile_flag=true
compile_flag=false
asan_flag=false
# run
run_flag=true
#run_flag=false
# plot
plot_flag=true
show_flag=true
#show_flag=false

### python options
# setOutput
outxyz=0
savejustelementxyz=0
outxyz_fname=$dir_run/input_epairs.xyz

#####################################################################################################################

mkdir -p $dir_run

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

### run
if ( $run_flag ) ; then
    for ile in ${!Lepairss[@]} ; do
	lepairs=${Lepairss[$ile]}
	if [ ${vdw:0:2} == 'LJ' ] ; then
	    longname="vdW$vdw-Hcorr$hbond-Ep$epairs_flag-LEpairs$lepairs"
	else
	    longname="vdW$vdw-kMorse$kmorse-Hcorr$hbond-Ep$epairs_flag-LEpairs$lepairs"
	fi
	echo ----------------------------------------------------
	echo $longname
	echo ----------------------------------------------------
	
	rm -f $dir_run/$longname.dat
	for (( istep=0 ; istep <= nsteps ; istep++ )) ; do
	    echo "  step = $istep / $nsteps"
	    lambda=`echo "$istep / $nsteps" | bc -l`
	    # write dofSelection.dat files
	    rm -f $dir_run/dofSelection.dat
	    touch $dir_run/dofSelection.dat
	    for idof in ${!dofnames[@]} ; do
		par=`echo "${dofs1[$idof]} + $lambda * ( ${dofs2[$idof]} - ${dofs1[$idof]} )" | bc -l`
		grep val_${dofnames[$idof]} $dir_data/dofSelection_run.dat | sed "s?val_${dofnames[$idof]}?$par?" >> $dir_run/dofSelection.dat
	    done
	    # run
	    python3 -u rigidscan.py \
		    --ivdw $ivdw --icoul $icoul --ihbond $ihbond --epairs $epairs --iepairs $iepairs --lepairs $lepairs --kmorse $kmorse \
		    --dof_selection $dir_run/dofSelection.dat --inputs_dir $dir_inputs --inputs $files_inputs \
		    --regularize $regularize --regcountweight $regcountweight \
		    --softclamp $softclamp --softclamp_start $softclamp_start --softclamp_max $softclamp_max \
		    --user_weights $user_weights --n_before $n_before --weight_a $weight_a --weight_alpha $weight_alpha --emin_min $emin_min --emin0 $emin0 \
		    --out_dir $dir_run --outxyz $outxyz --savejustelementxyz $savejustelementxyz --outxyz_fname $outxyz_fname 1> $dir_run/x.0
	    echo $lambda `tail -1 $dir_run/x.0` >> $dir_run/$longname.dat
	    rm $dir_run/dofSelection.dat $dir_run/x.0
	done
    done
fi	    

### plots
if ( $plot_flag ) ; then
    if [ ${vdw:0:2} == 'LJ' ] ; then
	longname="vdW$vdw-Hcorr$hbond-Ep$epairs_flag"
    else
	longname="vdW$vdw-kMorse$kmorse-Hcorr$hbond-Ep$epairs_flag"
    fi
    echo "set terminal pngcairo noenhanced crop size 1150,900 font 'Arial,18'" > $dir_run/x.gp
    echo "set output '$dir_run/$longname.png'" >> $dir_run/x.gp
    echo "set title '$longname'" >> $dir_run/x.gp
    echo "set xlabel 'Lambda'" >> $dir_run/x.gp
    echo "set ylabel 'Penalty function (eV^2)'" >> $dir_run/x.gp
    echo "plot "'\' >> $dir_run/x.gp
    for (( ile=0 ; ile < nLepairss-1 ; ile++ )) ; do
	echo "'$dir_run/$longname-LEpairs${Lepairss[$ile]}.dat'         u 1:2 w l lw 3 dt 1 lc $((ipar+1)) title 'Lepairs = ${Lepairss[$ile]} Ang',"'\' >> $dir_run/x.gp
    done
    echo "    '$dir_run/$longname-LEpairs${Lepairss[$nLepairss-1]}.dat' u 1:2 w l lw 3 dt 1 lc $nLepairss  title 'Lepairs = ${Lepairss[$nLepairss-1]} Ang'" >> $dir_run/x.gp
    gnuplot $dir_run/x.gp
    rm $dir_run/x.gp
    if ( $show_flag ) ; then eog $dir_run/$longname.png & fi
fi

exit
