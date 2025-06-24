#pragma once

//#ifndef  globals_h

//extern 
static int  verbosity = 1;
static int  idebug    = 0;

constexpr static const int ntmpstr=1024;
static char tmpstr[ntmpstr];

static double tick2second=1e-9;

// depending on debug / optimization leval 
#ifdef DEBUGBUILD
#define _assert( cond, action ) \
    if( !(cond) ){ \
        printf("Assertion failed: %s, line: %d  function: %s file: %s \n", #cond, __LINE__, __FUNCTION__, __FILE__ ); \
        printf("  => execute action: %s", #action ); \
        {action;} \
    }
#else
#define _assert( cond, action )
#endif


//#endif
