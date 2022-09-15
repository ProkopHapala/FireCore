
IF(WIN32)
    FIND_PATH( 
        OPENCL_SDK 
        NAMES CL/cl.h OpenCL/common/inc/CL/cl.h
        DOC "Path to the OpenCL SDK (Nvidia GPU Computing / AMD APP SDK)"
    )
    IF(NOT EXISTS ${OPENCL_SDK})
        MESSAGE(STATUS "OpenCL SDK (Nvidia GPU Computing / AMD APP SDK) not found! Please set Option 'OPENCL_SDK' if you need OpenCL!")
    ENDIF()
    
    FIND_PATH(
        OPENCL_INCLUDE_DIR
        NAMES CL/cl.h
        PATHS ${OPENCL_SDK}
        PATH_SUFFIXES include OpenCL/common/inc
    )
    
    IF(CAMPVIS_WIN32)
        FIND_LIBRARY(
            OPENCL_LIBRARY
            NAMES OpenCL.lib
            PATHS ${OPENCL_SDK}
            PATH_SUFFIXES lib/Win32 OpenCL/common/lib/Win32
        )
    ELSEIF(CAMPVIS_WIN64)
        FIND_LIBRARY(
            OPENCL_LIBRARY
            NAMES OpenCL.lib
            PATHS ${OPENCL_SDK}
            PATH_SUFFIXES lib/x64 OpenCL/common/lib/x64
        )
    ENDIF()

ELSEIF(UNIX)

    FIND_PATH(
        OPENCL_INCLUDE_DIR
        NAMES CL/cl.h
        PATHS ${OPENCL_INCLUDE_DIR}
              ${OPENCL_DIR}/include
              ${OPENCL_DIR}/OpenCL/common/inc
              $ENV{OPENCL_INCLUDE_DIR}
              $ENV{OPENCL_DIR}/include
              $ENV{OPENCL_DIR}/OpenCL/common/inc
              /usr/local/cuda/include
              /usr/local/include
              /usr/include
    )

    FIND_LIBRARY(
        OPENCL_LIBRARY
        NAMES OpenCL
        PATHS ${OPENCL_LIBRARY_DIR}
              ${OPENCL_DIR}/lib
              ${OPENCL_DIR}/lib/x86
              $ENV{OPENCL_LIBRARY_DIR}
              $ENV{OPENCL_DIR}/lib
              $ENV{OPENCL_DIR}/lib/x86
              /usr/local/lib64
              /usr/local/lib
              /usr/lib64
              /usr/lib
    )

ENDIF()

IF(OPENCL_INCLUDE_DIR AND OPENCL_LIBRARY)
    SET(OPENCL_FOUND TRUE)
    MARK_AS_ADVANCED(OPENCL_INCLUDE_DIR OPENCL_LIBRARY)
ELSE()
    SET(OPENCL_FOUND FALSE)
    MESSAGE(STATUS "OpenCL library not found (OPENCL_INCLUDE_DIR and/or OPENCL_LIBRARY missing)!")
    MARK_AS_ADVANCED(CLEAR OPENCL_INCLUDE_DIR OPENCL_LIBRARY)
ENDIF()
