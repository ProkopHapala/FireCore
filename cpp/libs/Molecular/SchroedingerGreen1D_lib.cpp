
#include <globals.h>

#include "SchroedingerGreen1D.h"

SchroedingerGreen1D W;

#include "libUtils.h"

extern "C"{

void setVerbosity( int verbosity_, int idebug_ ){
    verbosity = verbosity_;
    idebug    = idebug_;
}

void init_buffers(){
    // double*  V     = 0; // potential
    // double*  psi   = 0; // psi
    // double* fpsi   = 0; // derivatives of 
    // double* source = 0;
    buffers.insert( { "V",       W.V      } );
    buffers.insert( { "psi",     W.psi    } );
    buffers.insert( { "fpsi",    W.fpsi   } );
    buffers.insert( { "Apsi",    W.Apsi   } );
    buffers.insert( { "source",  W.source } );
    buffers.insert( { "EQF",    &W.E      } );
    
}

double  init( int nx, int ny, double dstep, double m_Me ){ return W.init(nx,dstep,m_Me); }
double stepResonance( double E0, double dt ){ return W.stepResonance( E0, dt ); };

} // extern "C"
