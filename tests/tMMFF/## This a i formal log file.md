##  This is a informal log file 
18 Feb

    Today all the expected Fix about p.add/sub is clrified and for the time being shift0 will be 0 and the p.sub(gridff.pos0) is final. I was making mistake due to the ignorance of vector notation.
    Now no shift no add only p.sub(gridff.pos0). 

    ##19 Feb
    The fix of symmetrization is pending and also its use.

    Need to chcek while calculating the morse or coulomb interaction which is taken apos_ (symmetrised) or apos (unsymmetrised), though there is two flag for this one is to make it symmetric
    which is called bSymterise, and the other one will be set to true, if the user wants to make it symmetric, which is called bSymmetrised. Whis is not mostly wotring.

    So today I can make the symmetrization function to work by either where it is defined set it to do nothing or use the apos while calculating all the interactions. 


    So today all the problem of single atom z scan comapred with LAMPS is solved.
    Points to be noted 
        1. The mismatch between LMAPS and FireCore zscan for single atom was not matching due to copule of things

        (a) the alpha value for the Morse pottential used in FireCore is 1.5 which was not trhe same as in LMAPS 
        (b) the there is no shift0 acting now in FireCore, the ztop of the substarte atom is same as the ztop of the moleculae atom then start the scan and will get the plot accordingly
        after some distance prefrebaly beyond just before R0. 
        
        2. The results are identitical in Firecore for both the codition with and without Grid, I mean the direct calculation with out forming the grid and the grid calculation.

        Now We need to check the xy scan for both Single atom and NaCl slab substrate by keeping the z at R0 from the z scan may be. Also need to figureit out about the exact posistion 
        of the molecule over slab.

    Today I am loging out 

    ## 20 Feb
    Today I will do the xy scan for single atom and NaCl slab substrate. Followed by PTCDA
    one more thing to check Grid and without Grid calculation for NaCl slab substrate. 
    Also eed to check about Columb interaction for the Grid calculation. 

    Today Only compared the 2d xy scan nothing more, this is realy very less work

    Need to implement the cutoff for the grid calculation. As it is there in LAMMPS.

    Also need to explore the Ewald sum for the grid calculation. how to make it work and all.

    mostly this columb is showing value 0 due to the +- 0.7 charge of NaCl substrate and no charge in C.iz0.xyz.
    
    I am signing out for the day.

    ## 21 Feb

    Today Lets try to do the solve the direct calculation without generating the grid and then compare it with the grid calculation. It was calling the symmetrization function, which has a problem.

    Yes it is done without symmetrization it is workign fine fised at GridFF.h:1250.

    Now lets think about the cutoff implemetation for the direct calculations. 


    In PLQ job it is not working fine with the Ewald sum callculation. the job mode Ewald is working instead which will save the umpy file with _ocl suffix. 


    ## 22 Feb 27 Feb

    I was debuging the PLQ job for the Ewald sum calculation. and comparing the with LAMMPS. WHICH IS NOW WORKING FINE. WITH THE OnLY CHAGE IS IN THE ORDER OF THE CALCULATION OF mORSE AND COULOMB potential.

    ## 28 Feb
    Today lets try the ocl in Firecore.

    ## 3 Mar


    ## 10 Mar
    Charge Analysis of NaCl slab substrate. using GPAW  
    Potential Fitting of NaCl slab substrate from GPAW.
    Relaxed scan of NaCl slab substrate.

    Tomorrow figure out the relaxed scan and required properties, 

    ## 11 Mar
    Did all the scan for NaCl molecule on top of NaCl slab substrate. 
    Need to do the same for PTCDA.

    ## 12 Mar
    Need to learn LAMMPS 
    
    Need to do theRelaxed scan of PTCDA on top of NaCl slab substrate.
    debuged the relaxed scan issue not giving the correct position of the molecule after realaxation. I was stupid did not trun on the MMFF flag as 1 it was set to -1 . Now that is working fine need to do the relaxed scan for PTCDA on top of NaCl slab substrate. 

    Need to do the scan for PTCDA on top of NaCl slab substrate.

    ## 13 Mar
    Need to do the rigid scan for PTCDA on top of NaCl slab substrate first and compare the results with LAMMPS, prepare final results and put it in the paper. 

    then may be the relaxed scan for PTCDA on top of NaCl slab substrate.

    ## 14 Mar
    Did the rigid scan for PTCDA on top of NaCl slab substrate.

    Updated the paper with results but the graphs are not production ready, need to do the final graphs, with proper line width, colour, legend text size and x and y labels. 


    Need to do the relaxed scan for PTCDA on top of NaCl slab substrate.


    ## 24 Mar
    Attending CI/CD workshop by IT4I 
    the constrained_scan is developed now need to do the uff part , which is implemented in ru_omp, need to add these scan constrain there in run_omp.













