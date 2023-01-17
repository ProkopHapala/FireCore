
#ifndef  SchroedingerGreen2D_h
#define  SchroedingerGreen2D_h


/*

E0 ... target energy
E  ... current Energy
L  ... is laplace operator

minimize{ (E-E0)^2 } by yi
(d/dyi) * (E-E0)^2    =     2*(E-E0) * dE/dyi 

E = <Y|L-V|Y> / <Y|Y>


L[ix  ,iy]*y[ix  ,iy] = y[ix  ,iy]* ( - 4*y[ix  ,iy]  + y[ix+1,iy] + y[ix-1,iy] + y[ix  ,iy+1] + y[ix  ,iy-1] )
L[ix+1,iy]*y[ix+1,iy] = y[ix+1,iy]* ( - 4*y[ix+1,iy]  + y[ix+2,iy] + y[ix  ,iy] + y[ix+1,iy+1] + y[ix+1,iy-1] )

dL[ix,iy]/dy[ix  ,iy] =  -8*y[ix,iy] +     (y[ix+1,iy] + y[ix-1,iy] + y[ix,iy+1] + y[ix,iy-1])
dL[ix,iy]/dy[ix+1,iy] =  +y[ix,iy]

dL[ix+1,iy  ]/dy[ix,iy] =  +y[ix+1,iy  ]
dL[ix-1,iy  ]/dy[ix,iy] =  +y[ix-1,iy  ]
dL[ix  ,iy+1]/dy[ix,iy] =  +y[ix  ,iy+1]
dL[ix  ,iy-1]/dy[ix,iy] =  +y[ix  ,iy-1]

dL/dy[ix,iy] = -8*y[ix,iy] +    2*(y[ix+1,iy] + y[ix-1,iy] + y[ix,iy+1] + y[ix,iy-1])

*/

class SchroedingerGreen2D{
    double dstep=1.;
    int nx=0,ny=0;
    double*  V  =0; // potential
    double*  psi=0; // psi
    double* fpsi=0; // derivatives of 

double sumE(){ // energy (hamiltonian)
    double E  =0;
    double rho=0;
    double invd2 = 1/(d*d);
    for(int ix=0; ix<nx; ix++){
        int i0y = nx+ix;
        for(int iy=0; iy<nx; iy++){
            int i = i0y+ix;
            double yi = psi[i];
            double vi = V  [i];
            //---- Laplacian
            // https://en.wikipedia.org/wiki/Discrete_Laplace_operator#Finite_differences
            double li = -4*yi;
            if((ix+1)<nx){ li+=psi[i+1 ]; }
            if((ix-1)>0 ){ li+=psi[i-1 ]; }
            if((iy+1)<ny){ li+=psi[i+nx]; }
            if((iy-1)>0 ){ li+=psi[i-nx]; }
            li*= invd2;
            // ---- Total energy
            E   += yi *( vi*yi* + li );
            rho += yi*yi;
        }
    }
}

double sumF( double renorm=-1 ){ // derivatives of energy
    double invd2 = 1/(dstep*dstep);
    for(int ix=0; ix<nx; ix++){
        int i0y = nx+ix;
        for(int iy=0; iy<nx; iy++){
            int i = i0y+ix;
            double yi = psi[i];
            double vi = V  [i];
            //---- Laplacian
            // https://en.wikipedia.org/wiki/Discrete_Laplace_operator#Finite_differences
            double dli = -4*yi;
            if((ix+1)<nx){ dli+=psi[i+1 ]; }
            if((ix-1)>0 ){ dli+=psi[i-1 ]; }
            if((iy+1)<ny){ dli+=psi[i+nx]; }
            if((iy-1)>0 ){ dli+=psi[i-nx]; }
            dli*= 2*invd2;
            // ---- Total energy
            fpsi[i] = (dli + 2*yi*vi)*renorm;
        }
    }
}


};

#endif

