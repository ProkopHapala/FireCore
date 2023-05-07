find_package( OpenCL REQUIRED )
#find_package( CL     REQUIRED )
find_package( clFFT  REQUIRED )

include_directories( ${COMMON_SRCS}/OpenCL    )
#include_directories( ${COMMON_SRCS}/molecular )

add_library(OCL_GridFF  SHARED OCL_GridFF.cpp )
target_link_libraries(OCL_GridFF ${OPENCL_LIBRARY} ${CLFFT_LIBRARIES})

add_library(MMFFmulti_lib  SHARED MMFFmulti_lib.cpp $<TARGET_OBJECTS:DynamicOpt> )
target_link_libraries(MMFFmulti_lib ${OPENCL_LIBRARY} ${CLFFT_LIBRARIES})