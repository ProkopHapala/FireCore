#find_package( OpenCL REQUIRED )
find_package( CL     REQUIRED )
find_package( clFFT  REQUIRED )

include_directories( ${COMMON_SRCS}/OpenCL    )
#include_directories( ${COMMON_SRCS}/molecular )

add_library(OCL_GridFF  SHARED OCL_GridFF.cpp ${CL_LIBRARY} ${clFFT_LIBRARY} )