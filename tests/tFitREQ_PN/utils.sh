function print_output() {
    local dir_run=$1
    local longname=$2
    local ihbond=$3
    local iepairs=$4
    grep -B3 'FitREQ::printDOFvalues()' $dir_run/dat/run-$longname.out | head -3
    if [ $iepairs -eq 0 ] ; then
	grep -B3 'FitREQ::printDOFvalues()' $dir_run/dat/run-$longname.out | tail -2
    else
	if [ $ihbond -eq 0 ] ; then
	    grep -A5 'FitREQ::printDOFvalues()' $dir_run/dat/run-$longname.out | tail -4
	else
	    grep -A7 'FitREQ::printDOFvalues()' $dir_run/dat/run-$longname.out | tail -6
	fi
    fi
    grep AFTER_REG $dir_run/dat/run-$longname.out | grep fDOF | awk '{print $4,$6,$8}' > $dir_run/dat/traj-$longname.dat
}

function plot_polarmap() {
    local dir_run=$1
    local out=$2
    local title=$3
    local dmin=$4
    local dmax=$5
    local Emin=$6
    local Emax=$7
    local f=$8
    local longname=$9
    echo "set terminal pngcairo enhanced crop size 500,900 font 'Arial,18'" > x.gp
    echo "set output '$dir_run/$out.png'" >> x.gp
    echo "set title '$title'" >> x.gp
    echo "set view map" >> x.gp
    echo "set pm3d interpolate 5,5" >> x.gp
    echo "set xlabel 'Angle (deg)'" >> x.gp
    echo "set xrange [-90:90]" >> x.gp
    echo "set xtics 30" >> x.gp
    echo "set ylabel 'Distance (Ang)'" >> x.gp
    echo "set yrange [$dmin:$dmax]" >> x.gp
    #echo "set logscale y" >> x.gp
    #echo "set ytics 1 nologscale" >> x.gp
    echo "set cblabel 'Energy (kcal/mol)'" >> x.gp
    echo "set cbrange [$Emin:$Emax]" >> x.gp
    echo "set palette defined (-1 'blue', 0 'white', 1 'red')" >> x.gp
    #echo "set size 0.95,1" >> x.gp
    echo "splot '$dir_run/dat/$f-$longname.dat' u 1:2:3 w pm3d notitle" >> x.gp
    gnuplot x.gp
    rm x.gp
}
