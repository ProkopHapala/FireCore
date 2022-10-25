
#ifndef MolWorld_sp3_ocl_h
#define MolWorld_sp3_ocl_h

#include "MolWorld_sp3.h"
#include "OCL_DFT.h"
#include "OCL_PP.h"

class MolWorld_sp3_ocl : public MolWorld_sp3 { public:
  OCL_PP     ocl;

// ======== Functions
void init_ocl(){
    printf( "DEBUG init_ocl() \n" );
    ocl.initPP( "common_resources/cl" );
    printf( "DEBUG init_ocl() END\n" );
}

void init( bool bGrid ){
    MolWorld_sp3::init(bGrid);
    init_ocl();
    printf( "DEBUG MolWorld::init() DONE \n");
}

};

#endif
