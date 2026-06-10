#!/bin/bash

### sysnames
#sysnames='C4H3NO2-C4H3NO2'
#sysnames='CH2NH-CH2NH'
#sysnames='H2O-H2O'
#sysnames='HBr-HBr'
#sysnames='HCl-HCl'
#sysnames='HCN-HCN'
#sysnames='HCONH2-HCONH2'
#sysnames='HCOOH-HCOOH'
#sysnames='HF-HF'
#sysnames='NH3-NH3'
sysnames='H2O-A1_H2O-D1-y'
#sysnames='CH2O-NH3'

### current working directory
cwdir=01-run

### kMorse
#kMorses='1.6 1.7 1.8 -1' # value of alpha for Morse and Buck interactions, a negative value means that alpha = 6 / R0
kMorses='1.7'

### setModel
#VDWs=( 'LJ' 'LJr8' 'LJr9' 'Morse' 'Buck') 
#iVDWs=('1'  '2'    '3'    '4'     '5'   )
VDWs=( 'Morse') 
iVDWs=('4'    )
iCOUL=1
#HBs=( 'nil' 'H1' 'H2')
#iHBs=('0'   '1'  '2' )
HBs=( 'H2')
iHBs=('2')
#EPs=( 'nil' 'SR' 'SR2' 'SR3' 'SR4')
#iEPs=('0'   '1'  '2'   '3'   '4'  )
EPs=( 'SR4')
iEPs=('4'  )

### LEpairs
#Lepairss='0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 2.0'  # distance of the Epair from the host atom
Lepairss='0.5'

### SR4
#SR4cuts='1.0 1.5 2.0' # cutoff for SR4 Epair interactions
SR4cuts='1.5'
SR4n=2
SR4m=2

### setPenalty
clamp=1 # Hardly restrain the values of parameters during optimization, 0=no, 1=yes
regularize=0
regcountweight=0
softclamp=1
softclamp_start=0.0867 # eV ~ 2 kcal; 0.0867 eV ~ 2 kcal, 0.1735 eV ~ 4 kcal, 0.2602 eV ~ 6 kcal
softclamp_max=0.1735 # eV ~ 4 kcal
user_weights=1
weight_a=1.0
weight_alpha=4.0

### setOptimization
ialg='2' # 0=move_GD, 1=move_MD, 2=move_GD_BB_short, 3=move_GD_BB_long
nstep=2000
fmax=1e-8

### folders
dir_main=/home/niko/work/HBOND/FireCore/tests/tFitREQ_PN
dir_firecore=/home/niko/work/HBOND/FireCore
dir_inputs=/home/niko/work/HBOND/REFERENCE/2-pairs_small_small/4-to_firecore/confs_wb97m
dir_data=$dir_main/data

#####################################################################################################################

icpu=1

source $dir_main/01-set_inputs.sh
source $dir_main/01-init_dofs.sh
source $dir_main/01-iterate_dofs.sh

### include dir_firecore into PYTHONPATH
export PYTHONPATH="$dir_firecore:$PYTHONPATH"

### compile Firecore
cd $dir_firecore/cpp/Build
make FitREQ_PN_lib || { echo "ERROR: FireCore compilation failed!"; exit 1; }
cd $dir_main

### loop over systems
for sysname in $sysnames ; do 

    ## inputs
    set_inputs $sysname
    echo "sysname=$sysname Files to be used for the fitting: $files_inputs"

    ## count runs
    nruns=0
    for iv in ${!VDWs[@]} ; do
	iVDW=${iVDWs[$iv]}
	kMorses_tmp=$kMorses ; if [ $iVDW -lt 4 ] ; then kMorses_tmp=1.8 ; fi
	for kMorse in $kMorses_tmp ; do
	    for ih in ${!HBs[@]} ; do
		iHB=${iHBs[$ih]}
		for ie in ${!EPs[@]} ; do
		    iEP=${iEPs[$ie]}
		    Lepairss_tmp=$Lepairss ; if [ $iEP -eq 0 ] ; then Lepairss_tmp=0.5 ; fi
		    for Lepairs in $Lepairss_tmp ; do
			SR4cuts_tmp=$SR4cuts ; if [ $iEP -ne 4 ] ; then SR4cuts_tmp=1.0 ; fi
			for SR4cut in $SR4cuts_tmp ; do
			    ((nruns++))
			done
		    done
		done
	    done
	done
    done
    echo "sysname=$sysname Number of different models to be fitted: $nruns"

    ## main loop
    irun=0
    for iv in ${!VDWs[@]} ; do
	VDW=${VDWs[$iv]}
	iVDW=${iVDWs[$iv]}
	kMorses_tmp=$kMorses ; if [ $iVDW -lt 4 ] ; then kMorses_tmp=1.8 ; fi
	for kMorse in $kMorses_tmp ; do
	    for ih in ${!HBs[@]} ; do
		HB=${HBs[$ih]}
		iHB=${iHBs[$ih]}
		for ie in ${!EPs[@]} ; do
		    EP=${EPs[$ie]}
		    iEP=${iEPs[$ie]}
		    epairs=$iEP ; if [ $iEP -gt 0 ] ; then epairs=1 ; fi
		    Lepairss_tmp=$Lepairss ; if [ $iEP -eq 0 ] ; then Lepairss_tmp=0.5 ; fi
		    for Lepairs in $Lepairss_tmp ; do
			SR4cuts_tmp=$SR4cuts ; if [ $iEP -ne 4 ] ; then SR4cuts_tmp=1.0 ; fi
			for SR4cut in $SR4cuts_tmp ; do
			    
			    longname="vdW$VDW-kMorse$kMorse-Hcorr$HB-Ep$EP-LEpairs$Lepairs-SR4cut$SR4cut"
			    irun=$(( irun + 1 ))
			    echo "-------------------------------------------------------------------------"
			    echo "sysname=$sysname $irun/$nruns $longname"
			    echo "-------------------------------------------------------------------------"
			    
			    # set initial values/ranges for fitting and count runs with different initial dof values
			    init_dofs $sysname $iHB $iEP
			    echo "sysname=$sysname Number of combinations of the initial values of the parameters: $dof_nruns"
			    if [ $dof_nruns -eq 0 ] ; then continue ; fi
			
			    # create working directory
			    dir_run=$dir_main/run_$sysname/$cwdir/runCPU.$icpu
			    mkdir -p $dir_run
			    cd $dir_run
			    
			    # run with all possible initial values
                            iterate_dofs

			    # clean-up working directory
			    rm -fr $dir_run

			done
		    done
		done
	    done
	done
    done
    
done

exit
