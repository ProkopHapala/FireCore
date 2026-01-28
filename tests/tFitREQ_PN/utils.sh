function count_runs() {
    nruntot=0
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
			nruntot=$(( nruntot + 1 ))
		    done
		done
	    done
	done
    done
}

function atom_types_begin_add() {
    local name=$1
    local unit=$2
    parnames=("${parnames[@]}" "$name")
    par=''
    for idof in ${!dofnames[@]} ; do
	if [ "${dofnames[$idof]}" == "$name" ] ; then
	    par=${dofstarts[$idof]}
	    break
	fi
    done
    if [ "$par" == '' ] ; then echo "ERROR $name NOT FOUND" ; exit; fi
    startpars=("${startpars[@]}" "$par")
    parunits=("${parunits[@]}" "$unit")
}

function atom_types_begin() {
    atom_types=()
    parnames=()
    startpars=()
    parunits=()
    # consider just atom types to be fitted
    atom_types_list=''
    for (( ia=0 ; ia < natom_types_all ; ia++ )) ; do
	at=${atom_types_all[$ia]}
	fit=${fits[$ia]}
	if ( ! $fit ) ; then continue ; fi
	atom_types=("${atom_types[@]}" "$at")
	atom_types_list="$atom_types_list $at"
	if [ $ihbond -gt 0 ] ; then
	    atom_types_begin_add "${at}.H" "-"
	fi
	if [ $iepairs -gt 0 ] ; then
	    atom_types_begin_add "E_${at}.R" "Ang"
	    atom_types_begin_add "E_${at}.H" "eV^1/2"
	fi
    done
    natom_types=${#atom_types[@]}
    npars=${#startpars[@]}
}    

function atom_types_end() {
    unset atom_types
    unset parnames
    unset startpars
    unset parunits
}
    
function autoinputs() {
    # find molecules
    mols_to_fit=''
    for mol in $mols ; do
	types=`head -2 $dir_mols/$mol.xyz | tail -1`
	for a in $atom_types_nofit ; do types=`echo $types | sed "s?$a??g"` ; done
	for at in "${atom_types[@]}" ; do
	    types=`echo $types | sed "s?$at??g"`
	done
	if [ "${types// }" == '' ] ; then mols_to_fit="$mols_to_fit $mol" ; fi
    done
    # find pairs
    files_inputs=''
    for pair_xyz in `ls $dir_inputs/*.xyz | sed "s?$dir_inputs/??g"` ; do
	pair=`echo $pair_xyz | sed "s?-? ?g" | sed "s?_? ?g" | awk '{print $1,$3}'`
	for mol in $mols_to_fit ; do
	    pair=`echo $pair | sed "s?$mol??g"`
	done
	if [ "${pair// }" == '' ] ; then files_inputs="$files_inputs $pair_xyz" ; fi
    done
}

function print_output() {
    grep -B3 'FitREQ_PN::printDOFvalues()' $dir_run/dat/run-$longname.out | head -3
    if [ $iepairs -eq 0 ] ; then
	if [ $ihbond -gt 0 ] ; then
	    grep -B3 'FitREQ_PN::printDOFvalues()' $dir_run/dat/run-$longname.out | tail -2
	fi
    else
	if [ $ihbond -eq 0 ] ; then
	    grep -A5 'FitREQ_PN::printDOFvalues()' $dir_run/dat/run-$longname.out | tail -4
	else
	    grep -A7 'FitREQ_PN::printDOFvalues()' $dir_run/dat/run-$longname.out | tail -6
	fi
    fi
    grep AFTER_REG $dir_run/dat/run-$longname.out | grep fDOF | awk '{print $4,$6,$8}' > $dir_run/dat/traj-$longname.dat
}

function store_output() {
    if [ ! -f $dir_run/${parnames[0]}.dat ] ; then
	echo "  NO $dir_run/${parnames[0]}.dat FILE PRESENT!"
	echo "  SOMETHING MUST HAVE GONE WRONG!"
	exit
    fi
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
}

function plot_traj() {
    echo "set terminal pngcairo noenhanced crop size 1150,900 font 'Arial,18'" > $dir_run/x.gp
    echo "set output '$dir_run/pars.png'" >> $dir_run/x.gp
    echo "set title 'Optimization trajectory'" >> $dir_run/x.gp
    echo "set xlabel 'Steps'" >> $dir_run/x.gp
    echo "set ylabel 'Parameter value'" >> $dir_run/x.gp
    echo "plot "'\' >> $dir_run/x.gp
    for (( ipar=0 ; ipar < npars-1 ; ipar++ )) ; do
	echo "'$dir_run/dat/${parnames[$ipar]}-$longname.dat'    u 1:2 w l lw 3 dt 1 lc $((ipar+1))  title '${parnames[$ipar]}',"'\' >> $dir_run/x.gp
    done
    echo "    '$dir_run/dat/${parnames[$npars-1]}-$longname.dat' u 1:2 w l lw 3 dt 1 lc $npars  title '${parnames[$npars-1]}'" >> $dir_run/x.gp
    gnuplot $dir_run/x.gp
    # traj
    echo "set terminal pngcairo noenhanced crop size 1400,900 font 'Arial,18'" > $dir_run/x.gp
    echo "set output '$dir_run/traj.png'" >> $dir_run/x.gp
    echo "set xlabel 'Steps'" >> $dir_run/x.gp
    echo "set ylabel 'Fitness (eV^2)' tc rgb 'black'" >> $dir_run/x.gp
    echo "set ytics nomirror tc rgb 'black'" >> $dir_run/x.gp
    echo "set logscale y" >> $dir_run/x.gp
    echo "set y2label 'Force norm (eV^2/a.u.^2)' tc rgb 'red'" >> $dir_run/x.gp
    echo "set y2tics nomirror tc rgb 'red'" >> $dir_run/x.gp
    echo "set logscale y2" >> $dir_run/x.gp
    echo "plot '$dir_run/dat/traj-$longname.dat' u 1:3 w l lw 3 lc rgb 'red'   axes x1y2 notitle,"'\' >> $dir_run/x.gp
    echo "     '$dir_run/dat/traj-$longname.dat' u 1:2 w l lw 3 lc rgb 'black' axes x1y1 notitle" >> $dir_run/x.gp
    gnuplot $dir_run/x.gp
    convert -background white -gravity North $dir_run/pars.png $dir_run/traj.png -append $dir_run/png/opt-$longname.png
    rm $dir_run/x.gp $dir_run/pars.png $dir_run/traj.png
}

function plot_cartmap() {
    local out=$1
    local title=$2
    local data=$3
    python3 -u plot_cartmap.py --filename $dir_run/dat/$data-$longname.dat --Emin $Emin --Emax $Emax --title "$title" --savefile $dir_run/cart.png --dmax $dmax
    convert $dir_run/cart.png -trim -border 10x10 -bordercolor white $dir_run/cart_$out.png
    rm $dir_run/cart.png
}    

function plot_polarmap() {
    local out=$1
    local title=$2
    local data=$3
    echo "set terminal pngcairo noenhanced crop size 500,900 font 'Arial,18'" > $dir_run/x.gp
    echo "set output '$dir_run/polar_$out.png'" >> $dir_run/x.gp
    echo "set title '$title'" >> $dir_run/x.gp
    echo "set view map" >> $dir_run/x.gp
    echo "set pm3d interpolate 5,5" >> $dir_run/x.gp
    echo "set xlabel 'Angle (deg)'" >> $dir_run/x.gp
    echo "set xrange [-90:90]" >> $dir_run/x.gp
    echo "set xtics 30" >> $dir_run/x.gp
    echo "set ylabel 'Distance (Ang)'" >> $dir_run/x.gp
    echo "set yrange [$dmin:$dmax]" >> $dir_run/x.gp
    #echo "set logscale y" >> $dir_run/x.gp
    #echo "set ytics 1 nologscale" >> $dir_run/x.gp
    echo "set cblabel 'Energy (kcal/mol)'" >> $dir_run/x.gp
    echo "set cbrange [$Emin:$Emax]" >> $dir_run/x.gp
    echo "set palette defined (-1 'blue', 0 'white', 1 'red')" >> $dir_run/x.gp
    echo "set size 0.95,1" >> $dir_run/x.gp
    echo "splot '$dir_run/dat/$data-$longname.dat' u 1:2:3 w pm3d notitle" >> $dir_run/x.gp
    gnuplot $dir_run/x.gp
    rm $dir_run/x.gp
}

function plot_minlines() {
    # rmin
    echo "set terminal pngcairo noenhanced crop size 700,900 font 'Arial,18'" > $dir_run/x.gp
    echo "set output '$dir_run/rmin.png'" >> $dir_run/x.gp
    #echo "set title '$f'" >> $dir_run/x.gp
    echo "set xlabel 'Angle (deg)'" >> $dir_run/x.gp
    echo "set xrange [-90:90]" >> $dir_run/x.gp
    echo "set xtics 30" >> $dir_run/x.gp
    echo "set ylabel 'Distance of the minimum (Ang)'" >> $dir_run/x.gp
    echo "set yrange [$dmin:$dmax]" >> $dir_run/x.gp
    echo "set logscale y" >> $dir_run/x.gp
    echo "set ytics 1 nologscale" >> $dir_run/x.gp
    echo "plot '$dir_run/dat/${f}__ref_lines-$longname.dat'   u 1:2 w lp lw 3 pt 6 ps 2 lc rgb 'black' title 'ref',"'\' >> $dir_run/x.gp
    echo "     '$dir_run/dat/${f}__model_lines-$longname.dat' u 1:2 w lp lw 3 pt 6 ps 2 lc rgb 'red'   title 'model'" >> $dir_run/x.gp
    gnuplot $dir_run/x.gp
    # Emin
    echo "set terminal pngcairo noenhanced crop size 700,900 font 'Arial,18'" > $dir_run/x.gp
    echo "set output '$dir_run/emin.png'" >> $dir_run/x.gp
    #echo "set title '$f'" >> $dir_run/x.gp
    echo "set xlabel 'Angle (deg)'" >> $dir_run/x.gp
    echo "set xrange [-90:90]" >> $dir_run/x.gp
    echo "set xtics 30" >> $dir_run/x.gp
    echo "set ylabel 'Energy of the minimum (kcal/mol)'" >> $dir_run/x.gp
    echo "plot '$dir_run/dat/${f}__ref_lines-$longname.dat'   u 1:3 w lp lw 3 pt 6 ps 2 lc rgb 'black' title 'ref',"'\' >> $dir_run/x.gp
    echo "     '$dir_run/dat/${f}__model_lines-$longname.dat' u 1:3 w lp lw 3 pt 6 ps 2 lc rgb 'red'   title 'model'" >> $dir_run/x.gp
    gnuplot $dir_run/x.gp
    convert -background white -gravity North $dir_run/rmin.png $dir_run/emin.png +append $dir_run/png/min_$f-$longname.png
    rm $dir_run/x.gp $dir_run/rmin.png $dir_run/emin.png
}

function find_eminmax() {
    local file=$1
    local de=$2
    grep -v nan $file | grep -v angle > $dir_run/x.0
    cd $dir_run
    Emin=`echo "x.0 $dmin $dmax" | $dir_run/minmax.x | awk '{print $1}'`
    rm x.0
    cd $dir_exec
    if (( $(echo "$de > 0.0" | bc -l) )) ; then
	Emax=$(awk -v emin="$Emin" -v de="$de" 'BEGIN{print emin+de}')
    else
	Emax=$(awk -v emin="$Emin" 'BEGIN{print -(emin)}')
    fi
    #awk -v emin="$Emin" -v emax="$Emax" 'BEGIN{if (!(emin<0 && emax>0)) exit 1}' || { echo "ERROR: Emin must be <0 and Emax >0 (got Emin=$Emin, Emax=$Emax)"; exit 1; }
}    

function find_eminmaxdiff() {
    local file=$1
    grep -v nan $file | grep -v angle > $dir_run/x.0
    cd $dir_run
    Emin=`echo "x.0 $dmin $dmax" | $dir_run/minmaxdiff.x | awk '{print $1}'`
    rm x.0
    cd $dir_exec
    Emax=$(awk -v emin="$Emin" 'BEGIN{print -(emin)}')
    #awk -v emin="$Emin" -v emax="$Emax" 'BEGIN{if (!(emin<0 && emax>0)) exit 1}' || { echo "ERROR: Emin must be <0 and Emax >0 (got Emin=$Emin, Emax=$Emax)"; exit 1; }
}    

function calc_rmse() {
    local file=$1
    grep -v nan $file | grep -v angle > $dir_run/x.0
    cd $dir_run
    rmse=`echo "x.0" | $dir_run/rmse.x`
    rm x.0
    cd $dir_exec
}    

function check_tbd_run() {
    local file=$1
    if [ -f $file ] ; then
	local out=`grep "save_grid_gnuplot(): saving to" $file`
    else
	local out=''
    fi
    if [ "$out" == '' ] ; then
	tbd=true
    else
	tbd=false
    fi
}

function check_tbd_plot() {
    local file=$1
    if [ -s $file ] ; then
	tbd=false
    else
	tbd=true
    fi
}
