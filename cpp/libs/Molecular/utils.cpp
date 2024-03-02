
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <vector>
#include <math.h>



extern "C"{



void  find_Alphabet( int n, double* Emap, double Etot_min, double Etot_max, double  minContrast ){


    // Find all alphabest of size 2
    int n=0;
    int level2[n*n];
    for( int i=0; i<n; i++ ){
        for( int j=0; j<i; j++ ){
            double Eij=Emap[i*n+j];
            double Eii=Emap[i*n+i];
            double Ejj=Emap[j*n+j];
            double contrast = Eij - (Eii+Ejj)*0.5;
            if( (Eij>Etot_min) && (Eij<Etot_max) && (  contrast>minContrast) ){
                level2[n  ]=i;
                level2[n+1]=j;
                int n++;
            }
        }
    }

    // find larger alphabets

    for( int i=0; i<n; i++ ){ // for all alphabets
        for( int j=0; j<i; j++ ){ // for all pairs (not in alphabet)


            for( int j=0; j<i; j++ ){ // for all pairs in alphabet

            }

            double Eij=Emap[i*n+j];
            double Eii=Emap[i*n+i];
            double Ejj=Emap[j*n+j];
            double contrast = Eij - (Eii+Ejj)*0.5;
            if( (Eij>Etot_min) && (Eij<Etot_max) && (  contrast>minContrast) ){
                level2[n  ]=i;
                level2[n+1]=j;
                int n++;
            }
        }
    }

}





} // extern "C"
