

constexpr int ntmpstr=2048;
char tmpstr[ntmpstr];

#include <globals.h>
//int verbosity = 1;
//int idebug    = 0;
//double tick2second=1e-9;

#include "SchroedingerGreen2D.h"

SchroedingerGreen2D W;

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

double  init( int nx, int ny, double dstep, double m_Me ){ return W.init(nx,ny,dstep,m_Me); }

double setStep      ( double dstep, double m_Me){ return W.setStep( dstep, m_Me ); }
double step         ( double E0, double dt ){ return W.step( E0, dt); };
double step_Green   (                      ){ return W.step_CG();     };
double stepResonance( double E0, double dt ){ return W.stepResonance( E0, dt ); };

int solve_Green( int maxIters, double maxErr ){ return W.solve_CG( maxIters, maxErr ); };

//double sumE(){ W.sumE(); }
//double sumF( double renorm1 ){ W.sumF(); }

} // extern "C"
