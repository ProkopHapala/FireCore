#pragma once

//#ifndef  globals_h

//extern
static int  verbosity = 1;
static int  idebug    = 0;
static int  id_DBG    = 1;

constexpr static const int EXCL_MAX = 16;
constexpr static const int ntmpstr=1024;
static char tmpstr[ntmpstr];

static double tick2second=1e-9;

static const double const_eVA2_Nm = 16.02176634;

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

#define DEBUG   printf( "DEBUG #l %i %s \n",    __LINE__, __FUNCTION__ );
#define DEBUGF  printf( "DEBUG #l %i %s %s \n", __LINE__, __FUNCTION__, __FILE__ );

//void dbg(char* s){ printf("DEBUG (%s) \n", s); };

#define DBG(format,args...) { printf("DEBUG "); printf(format,## args); }

//#endif
