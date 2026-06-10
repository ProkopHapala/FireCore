# this function takes sysname, iHB and iEP and sets all arrays of parameters for fitting
init_dofs() {
    local mysysname=$1
    local myiHB=$2
    local myiEP=$3
    # define atoms involved
    if [ "$mysysname" == 'C4H3NO2-C4H3NO2' ] ; then
	local atoms=('H_N' 'N_3' 'O_2')
    elif [ "$mysysname" == 'CH2NH-CH2NH' ] ; then
	local atoms=('H_N' 'N_2')
    elif [ "$mysysname" == 'H2O-H2O' ] ; then 
	local atoms=('H_O' 'O_3')
    elif [ "$mysysname" == 'HBr-HBr' ] ; then
	local atoms=('H_Br' 'Br')
    elif [ "$mysysname" == 'HCl-HCl' ] ; then
	local atoms=('H_Cl' 'Cl')
    elif [ "$mysysname" == 'HCN-HCN' ] ; then
	local atoms=('H_C1' 'N_1')
    elif [ "$mysysname" == 'HCONH2-HCONH2' ] ; then
	local atoms=('H_N' 'N_3' 'O_2')
    elif [ "$mysysname" == 'HCOOH-HCOOH' ] ; then
	local atoms=('H_O' 'O_2' 'O_3')
    elif [ "$mysysname" == 'HF-HF' ] ; then
	local atoms=('H_F' 'F_')
    elif [ "$mysysname" == 'NH3-NH3' ] ; then
	local atoms=('H_N' 'N_3')
    elif [ "$mysysname" == 'H2O-A1_H2O-D1-y' ] ; then 
	local atoms=('H_O' 'O_3')
    elif [ "$mysysname" == 'CH2O-NH3' ] ; then
	local atoms=('H_N' 'N_3' 'O_2')
    else
        echo "ERROR: $mysysname not implemented!"
        exit 1
    fi
    # set arrays with parameters to be written in dofSelection.dat
    dof_typenames=()
    dof_components=()
    dof_names=()
    dof_mins=()        # hard limit for min			
    dof_maxs=()	       # hard limit for max			
    dof_xlos=()	       # kick-in of lower regularization		
    dof_xhis=()	       # kick-in of higher regularization	
    dof_xstarts=()     # starting values for dofs, note: - if SR4, E_*.R and E_*.H are actually energies in eV
                       #                                 - otherwise, E_*.R is a distance in Ang and E_*.H is an energy in eV
    dof_klos=()	       # force constant for lower regularization	
    dof_khis=()	       # force constant for higher regularization
    dof_k0s=()	       # force constant for center regularization
    dof_imasss=()      # inverse of dof mass (used for MD)
    for ia in ${!atoms[@]} ; do
	local at=${atoms[$ia]}
	if [ $myiHB -gt 0 ] ; then
	    dof_typenames+=("$at")
	    dof_components+=('3')
	    dof_names+=("$at.H")
	    if [ "${at:0:2}" == 'H_' ] ; then
		dof_mins+=('1e-9')
		dof_maxs+=('5.0')
		dof_xlos+=('0.1')
		dof_xhis+=('0.9')
		dof_xstarts+=('.2 2')
	    else
		dof_mins+=('-5.0')
		dof_maxs+=('-1e-9')
		dof_xlos+=('-0.9')
		dof_xhis+=('-0.1')
		dof_xstarts+=('-.2 -2')
	    fi
	    dof_klos+=('0.0')
	    dof_khis+=('0.0')
	    dof_k0s+=('0.0')
	    dof_imasss+=('1.0')
	fi
	if [ $myiEP -gt 0 ] ; then
	    dof_typenames+=("E_$at" "E_$at")
	    dof_components+=('0' '3')
	    dof_names+=("E_$at.R" "E_$at.H")
	    if [ $myiEP -eq 4 ] ; then
		if [ "${at:0:2}" == 'H_' ] ; then
		    dof_mins+=('1e-9' '1e-9')
		    dof_maxs+=('50.0' '50.0')
		    dof_xlos+=('0.1' '0.1' )
		    dof_xhis+=('10.0' '10.0')
		    dof_xstarts+=('.2 2' '.2 2')
		else
		    dof_mins+=('-50.0' '-50.0')
		    dof_maxs+=('-1e-9' '-1e-9')
		    dof_xlos+=('-10.0' '-10.0')
		    dof_xhis+=('-0.1' '-0.1')
		    dof_xstarts+=('-.2 -2' '-.2 -2')
		fi
	    else
		if [ "${at:0:2}" == 'H_' ] ; then
		    dof_mins+=('1e-9' '1e-9')
		    dof_maxs+=('10.0' '50.0')
		    dof_xlos+=('0.1' '0.1' )
		    dof_xhis+=('2.0' '10.0')
		    dof_xstarts+=('.2 1' '.2 2')
		else
		    dof_mins+=('1e-9' '-50.0' )
		    dof_maxs+=('10.0' '-1e-9' )
		    dof_xlos+=('0.1' '-10.0' )
		    dof_xhis+=('2.0' '-0.1' )
		    dof_xstarts+=('.2 1' '-.2 -2')
		fi
	    fi
	    dof_klos+=('0.0' '0.0')
	    dof_khis+=('0.0' '0.0')
	    dof_k0s+=('0.0' '0.0')
	    dof_imasss+=('1.0' '1.0')
	fi
    done
    # count the number of runs
    if [ ${#dof_typenames[@]} -eq 0 ] ; then
	dof_nruns=0
	echo "WARNING: nothing to fit!"
    else
	dof_nruns=1
	for id in ${!dof_typenames[@]} ; do
	    local i=0
	    for ix in ${dof_xstarts[$id]} ; do
		((i++))
	    done
	    dof_nruns=$(( dof_nruns * i ))
	done
    fi
}
