# plot file TI_plot.dat with gnuplot. In the file, there are 3 columns "lamda", "FreeEnergy" and "Reference", which are the x and 2 y values. Add legend from the header of the file.
# /bin/python3 /home/kocimil1/Documents/testing_space/Free_energy/Three_particle_problem/main.py


gnuplot -e "
    set terminal pngcairo size 1920,1080 enhanced font 'Verdana,10';
    set output 'results/TI_plot.png';
    set title 'TI plot';
    set xlabel 'lamda';
    set ylabel 'Energy [eV]';
    set xtics 0.2;
    set format x '%.1f';
    set mxtics 5;

    plot 'results/TI_plot.dat' using 1:2:3 with yerrorbars title 'Thermodynamic integration', \
         'results/TI_plot.dat' using 1:6 with lines title 'Reference';
"
# show the plot
#display results/TI_plot.png

