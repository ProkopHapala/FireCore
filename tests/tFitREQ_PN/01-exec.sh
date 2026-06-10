#!/bin/bash

### passed from the script
icpu=$1

### passed from the list
irun=$2
VDW=$3
iVDW=$4
kMorse=$5
HB=$6
iHB=$7
EP=$8
iEP=$9
Lepairs=${10}
SR4cut=${11}

### passed from the fixed file
sysname=${12}
iCOUL=${13}
SR4n=${14}
SR4m=${15}
clamp=${16}
regularize=${17}
regcountweight=${18}
softclamp_start=${19}
softclamp_max=${20}
softclamp=${21}
user_weights=${22}
weight_alpha=${23}
weight_a=${24}
ialg=${25}
nstep=${26}
fmax=${27}
dir_main=${28}
dir_firecore=${29}
dir_inputs=${30}
dir_data=${31}

### needed to run the fitting
epairs=$iEP ; if [ $iEP -gt 0 ] ; then epairs=1 ; fi

### needed to store output
longname="vdW$VDW-kMorse$kMorse-Hcorr$HB-Ep$EP-LEpairs$Lepairs-SR4cut$SR4cut"

#####################################################################################################################

source $dir_main/01-init_dofs.sh
source $dir_main/01-set_inputs.sh
source $dir_main/01-iterate_dofs.sh

### set initial values/ranges for fitting and count runs with different initial dof values
init_dofs $sysname $iHB $iEP
echo "sysname=$sysname Number of combinations of the initial values of the parameters: $dof_nruns"
if [ $dof_nruns -eq 0 ] ; then rm -f CPU.$icpu ; exit ; fi
    
### include dir_firecore into PYTHONPATH
export PYTHONPATH="$dir_firecore:$PYTHONPATH"

### inputs
set_inputs $sysname
echo "sysname=$sysname Files to be used for the fitting: $files_inputs"

## create working directory
mkdir -p runCPU.$icpu
cd runCPU.$icpu
			    
## run with all possible initial values
iterate_dofs
			    
## clean-up working directory
cd ..
rm -fr runCPU.$icpu CPU.$icpu

exit
