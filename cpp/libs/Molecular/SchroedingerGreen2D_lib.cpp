

constexpr int ntmpstr=2048;
char tmpstr[ntmpstr];

int verbosity = 1;
int idebug    = 0;
//double tick2second=1e-9;

#include "SchroedingerGreen2D.h"

SchroedingerGreen2D W;

#include "libUtils.h"

extern "C"{

void init_buffers(){
    // double*  V     = 0; // potential
    // double*  psi   = 0; // psi
    // double* fpsi   = 0; // derivatives of 
    // double* source = 0;
    buffers.insert( { "V",       W.V      } );
    buffers.insert( { "psi",     W.psi    } );
    buffers.insert( { "fpsi",    W.fpsi   } );
    buffers.insert( { "source",  W.source } );
    buffers.insert( { "EQF",    &W.E      } );
}

void   init( int nx, int ny       ){ W.init(nx,ny); }
double step      ( double E0, double dt ){ return W.step( E0, dt); };
double step_Green(                      ){ return W.step_CG();     };

//double sumE(){ W.sumE(); }
//double sumF( double renorm1 ){ W.sumF(); }

} // extern "C"
