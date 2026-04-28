#! /bin/bash
grep BEFOR_REG OUT | grep " DOFs" | cut -c 62- > DOFs.log
gnuplot -persist -e "set title 'DOFs'; plot for [col=1:7] 'DOFs.log' using 0:col with lines title sprintf('Column %d', col)"


grep BEFOR_REG OUT | grep "fDOFs" | cut -c 62- > fDOFs.log
gnuplot -persist -e "set title 'fDOFs'; plot for [col=1:7] 'fDOFs.log' using 0:col with lines title sprintf('Column %d', col)"