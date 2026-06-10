# this function performs the fitting of a model for all possible initial values of the parameters
# it requires:
# - all arrays initialized by init_dofs
# - 01-fit.py script in dir_main with all files needed
# - irun to save reference data
# it puts output in tar archives in the parent directory

iterate_dofs() {
    # --- Staging local tracking variables
    local n_dofs=${#dof_typenames[@]}
    local val_arrays=()
    local val_counts=()
    
    for (( d=0 ; d < n_dofs ; d++ )) ; do
        local items=(${dof_xstarts[$d]})
        val_counts+=(${#items[@]})
    done

    # --- Synchronized Cartesian loop tracker
    local run_counter=0
    local current_indices=()
    for (( d=0 ; d < n_dofs ; d++ )) ; do current_indices+=(0) ; done

    while true ; do
        ((run_counter++))
	counter=$(printf "%05d" $run_counter)

	# Define a current xstart array
	local current_xstarts=()
	for (( d=0 ; d < n_dofs ; d++ )) ; do
            local idx=${current_indices[$d]}
	    local items=(${dof_xstarts[$d]})
	    current_xstarts+=(${items[$idx]})
	done

	# Unambiguous name for the run
	local unambiguous=${dof_typenames[0]}${current_xstarts[0]}
	for (( d=1 ; d < n_dofs ; d++ )) ; do
	    unambiguous="$unambiguous-${dof_typenames[$d]}${current_xstarts[$d]}"
	done

	# Skip if already ran
	if tar -xOf ../$longname.log.tar --wildcards --no-anchored "$unambiguous.log*" 2>/dev/null | grep -q "save_grid_gnuplot(): saving to" ; then
	    continue
	fi

	# Write dofSelection.dat files
	rm -f dofSelection.dat
	for (( d=0 ; d < n_dofs ; d++ )) ; do
	    echo "${dof_typenames[$d]} ${dof_components[$d]} ${dof_mins[$d]} ${dof_maxs[$d]} ${dof_xlos[$d]} ${dof_xhis[$d]} ${dof_klos[$d]} ${dof_khis[$d]} ${dof_k0s[$d]}  ${current_xstarts[$d]} ${dof_imasss[$d]}" >> dofSelection.dat
        done

	# Run
	echo "run $counter $unambiguous" > $unambiguous.log
	python $dir_main/01-fit.py \
	       --ivdw $iVDW --icoul $iCOUL --ihbond $iHB --epairs $epairs --iepairs $iEP --lepairs $Lepairs --kmorse $kMorse \
	       --sr4cut $SR4cut --sr4m $SR4m --sr4n $SR4n \
	       --dof_selection dofSelection.dat --inputs_dir $dir_inputs --inputs $files_inputs \
	       --fetypes $dir_data/ElementTypes.dat --fatypes $dir_data/AtomTypes.dat \
	       --clamp $clamp --regularize $regularize --regcountweight $regcountweight \
	       --softclamp $softclamp --softclamp_start $softclamp_start --softclamp_max $softclamp_max \
	       --user_weights $user_weights --weight_a $weight_a --weight_alpha $weight_alpha \
	       --ialg $ialg --nstep $nstep --fmax $fmax --out_dir . #>> $unambiguous.log
exit	
	# Store output 
	local Err=`grep "VERY FINAL |E|" $unambiguous.log | sed "s?=? ?" | awk '{print $4}'`
	local current_xends=()
	read -a current_xends < <(grep "VERY FINAL |E|" $unambiguous.log | sed "s?.*DOFs=\s*??")
	echo "#parnames $counter ${dof_names[*]}"            >  $unambiguous.fit
	echo "#initpars $counter ${current_xstarts[*]}"      >> $unambiguous.fit
	echo "#finalpars $counter ${current_xends[*]}"       >> $unambiguous.fit
	printf "#penalty %05d %.20f\n" "$run_counter" "$Err" >> $unambiguous.fit

	# Printout
	if ! grep -q CONVERGED $unambiguous.log ; then # not converged
	    printf "sysname=%s CPU %-4d Run %-6d NOT CONVERGED  Err = %16.12f kcal/mol\n" "$sysname" "$icpu" "$run_counter" "$Err"
	    for (( d=0 ; d < n_dofs ; d++ )) ; do
		printf "%48s INIT = %8.4f   FINAL = %8.4f\n" "${dof_names[d]}" "${current_xstarts[$d]}" "${current_xends[$d]}"
	    done
	    mv $unambiguous.log $unambiguous.lognc
	else # successfully finished
	    printf "sysname=%s CPU %-4d Run %-6d Steps = %-6d Err = %16.12f kcal/mol %s\n" "$sysname" "$icpu" "$run_counter" "`grep "^step=" $unambiguous.log | tail -1 | awk '{print $2}'`" "$Err"
	    for (( d=0 ; d < n_dofs ; d++ )) ; do
		printf "%48s INIT = %8.4f   FINAL = %8.4f\n" "${dof_names[d]}" "${current_xstarts[$d]}" "${current_xends[$d]}"
	    done
	fi
	
	# Save ref data only once
	if [[ $irun -eq 1 && $run_counter -eq 1 ]] ; then
	    for f in ${files_inputs//.xyz} ; do
		for suf in ref ref_lines ; do
		    mv ${f}__$suf.dat ../${f}__$suf.dat
		done
	    done
	fi
	
	# Create or insert into archives
	local flag=rf ; if [ $run_counter -eq 1 ] ; then flag=cvf ; fi
	tar -$flag ../$longname.log.tar $unambiguous.log* 1> /dev/null
	tar -$flag ../$longname.fit.tar $unambiguous.fit 1> /dev/null
	for f in ${files_inputs//.xyz} ; do
	    for suf in model model_lines ; do
		mv ${f}__$suf.dat $unambiguous-${f}__$suf.dat
		tar -$flag ../$longname-${f}__$suf.dat.tar $unambiguous-${f}__$suf.dat 1> /dev/null
	    done
	done
	
	# Cleanup folder
	rm -f *.dat $unambiguous.log* $unambiguous.fit $unambiguous-*__*.dat

        # --- Advance Cartesian tracker wheel odometer
        local dim=$(( n_dofs - 1 ))
        while [ $dim -ge 0 ] ; do
            ((current_indices[dim]++))
            if [ ${current_indices[$dim]} -lt ${val_counts[$dim]} ] ; then
                break
            fi
            current_indices[$dim]=0
            ((dim--))
        done

	# --- Exit condition: all parameter values done
        if [ $dim -lt 0 ] ; then
	    break
	fi
    done
}
