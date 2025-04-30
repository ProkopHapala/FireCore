#cmake_minimum_required ( VERSION 2.8 )
cmake_minimum_required(VERSION 3.10)

project ( SimpleSimulationEngine )

option(WITH_SDL    "build with GUI and 3D graphics (SDL+OpenGL)"                           OFF)
option(WITH_LUA    "build with Lua 5.2 for scripting inside interactive GUI"               OFF)
option(WITH_FFTW   "build with FFTW3 for scripting inside interactive GUI"                 OFF)
option(WITH_OPENCL "build with OpenCL GPU acceleration? up to 100x speedup"                OFF)
option(WITH_CUDA   "build with CUDA GPU acceleration?"                                     OFF)
option(WITH_OMP    "build with OpenMP paralelization for multiple CPUs"                    OFF)
option(WITH_AVX    "build with AVX256 SIMD intrinstic ? up to 10x speedup on 1CPU"         OFF)
option(RELEASE     "build optimized release (-O3/-Ofast) rather than debug? (-g -Og)"      OFF)
option(WITH_ASAN   "use runtime memory sanitizer (asan) to trace segfaults and mem-leaks?" OFF)

message("OPTIONS: -DWITH_SDL=${WITH_SDL} -DWITH_LUA=${WITH_LUA} -DWITH_FFTW=${WITH_FFTW} -DWITH_OPENCL=${WITH_OPENCL}  -DWITH_CUDA=${WITH_CUDA}  -DWITH_OMP=${WITH_OMP} -DWITH_AVX=${WITH_AVX} -DRELEASE=${RELEASE} -DWITH_ASAN=${WITH_ASAN}" )

# Enable CUDA language if selected
if(WITH_CUDA)
    # Use a supported GCC for nvcc host compilation to avoid version mismatch
    set(CMAKE_CUDA_HOST_COMPILER "/usr/bin/gcc-12")
    # Override nvcc version check for both identification and subsequent compilations
    set(CMAKE_CUDA_COMPILER_ID_FLAGS_INIT "${CMAKE_CUDA_COMPILER_ID_FLAGS_INIT} -allow-unsupported-compiler")
    set(CMAKE_CUDA_FLAGS_INIT "${CMAKE_CUDA_FLAGS_INIT} -allow-unsupported-compiler")
    project(SimpleSimulationEngine LANGUAGES C CXX) # Keep project simple first
    enable_language(CUDA)                          # Enable CUDA language support
    message(STATUS "CUDA language enabled.")
else()
    project(SimpleSimulationEngine LANGUAGES C CXX)
endif()

# ==============================
# ==== FLAGS
# ==============================

if (WITH_AVX)
    add_definitions(-DWITH_AVX)
endif()

if(WITH_LUA)
    add_definitions(-DWITH_LUA)
    find_package(Lua 5.2 REQUIRED)
    # Add more Lua-specific setup here, e.g., including directories, linking libraries
else()
    set(LUA_LIBRARIES "")
endif()

if(WITH_FFTW)
    add_definitions(-DWITH_FFTW)
    #find_package(FFTW3 REQUIRED)
    find_package(PkgConfig REQUIRED)
    pkg_search_module(FFTW3 REQUIRED fftw3)
    include_directories(${FFTW3_INCLUDE_DIRS})
    link_directories(${FFTW3_LIBRARY_DIRS})
    set(FFTW_LIBRARIES -lfftw3 -lfftw3_omp ${OpenMP_CXX_LIBRARIES})
    #set(FFTW_LIBRARIES -lfftw3 -lfftw3_omp )
else()
    set(FFTW_LIBRARIES "")
endif()


if(RELEASE)
    set(CMAKE_BUILD_TYPE Release)
else()
    set(CMAKE_BUILD_TYPE Debug)
endif()
if( UNIX )
    #SET( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -std=c++20" )
    SET( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++20 -fPIC -fno-strict-aliasing" )
endif()
# https://stackoverflow.com/questions/7724569/debug-vs-release-in-cmake
message ( "CMAKE_BUILD_TYPE "   ${CMAKE_BUILD_TYPE} )
if (  CMAKE_BUILD_TYPE MATCHES Debug )
    if ( WITH_ASAN )
        SET( AXULIARY_COMPILE_FLAGS "-g -Og -march=native -fsanitize=address -fno-omit-frame-pointer")
        #SET( AXULIARY_COMPILE_FLAGS "-g -Og -fsanitize=address -fno-omit-frame-pointer -static-libasan") # static-asan does not work see https://github.com/google/sanitizers/issues/796
    else()
        #SET( AXULIARY_COMPILE_FLAGS "-Wall -g -Og")
        ### flags for debug
        SET( AXULIARY_COMPILE_FLAGS "-g -Og -march=native")
    endif()
    #set(CMAKE_CXX_CPPCHECK "cppcheck")
    #message("DEBUG BUILD FFLAGS : " )
else()
    #SET( AXULIARY_COMPILE_FLAGS "-w -O3 -fPIC ")
    #SET( AXULIARY_COMPILE_FLAGS "-w -O2 -march=native")
    #SET( AXULIARY_COMPILE_FLAGS "-w -O2 -march=native")
    SET( AXULIARY_COMPILE_FLAGS "-w -Ofast -march=native -mtune=native")
    #SET( AXULIARY_COMPILE_FLAGS "-w -Ofast -march=native -fopt-info-vec-all")
endif()
if ( WITH_OMP )
    find_package(OpenMP)
    #if(OpenMP_CXX_FOUND)
    #    target_link_libraries(MyTarget PUBLIC OpenMP::OpenMP_CXX)
    #endif()
    SET( AXULIARY_COMPILE_FLAGS "${AXULIARY_COMPILE_FLAGS} -fopenmp" )
endif()

### flags always there
SET( IGNORE_WARRNING_FLAGS "-Wno-int-to-pointer-cast -Wno-switch -Wno-reorder -Wno-deprecated -Wno-maybe-uninitialized -Wno-char-subscripts -Wno-write-strings -Wno-format -Wno-parentheses -Wno-unused-but-set-variable  -Wno-narrowing -Wno-unused-result  -Wno-sign-compare -Wno-strict-aliasing -Wno-unused-variable -Wno-unused-value -Wno-comment -Wno-misleading-indentation " )
SET( WARRNING_TO_ERROR     "-Werror=return-type -Werror=init-self -Werror=uninitialized -Werror=implicit-fallthrough -Werror=tautological-compare -Werror=delete-non-virtual-dtor")
SET( CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${AXULIARY_COMPILE_FLAGS} ${IGNORE_WARRNING_FLAGS} ${WARRNING_TO_ERROR}" )
message("BUILD FFLAGS : " )
message ( ${CMAKE_CXX_FLAGS} )


# <<< CUDA Change Start >>>
if(WITH_CUDA)
    message(STATUS "Setting up CUDA flags...")
    # Set CUDA C++ Standard (must match CMAKE_CXX_STANDARD if device code uses C++ features)
    set(CMAKE_CUDA_STANDARD 20)
    set(CMAKE_CUDA_STANDARD_REQUIRED ON)
    set(CMAKE_CUDA_EXTENSIONS OFF) # Recommended

    # --- TODO CUDA: Set Target GPU Architectures ---
    # Option 1: Use CMAKE_CUDA_ARCHITECTURES (CMake 3.18+) - semicolon-separated list of compute capabilities (without 'sm_')
    # set(CMAKE_CUDA_ARCHITECTURES "70;75;86") # Example: Volta, Turing, Ampere (virtual arch only)
    # set(CMAKE_CUDA_ARCHITECTURES "70-real;75-real;86-real") # Example: Volta, Turing, Ampere (real arch only)
    # set(CMAKE_CUDA_ARCHITECTURES "75") # Example: Target only Turing

    # Option 2: Add flags directly (more compatible with older CMake)
    # Common practice: build for a range of architectures
    # set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -gencode arch=compute_70,code=sm_70") # Volta
    # set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -gencode arch=compute_75,code=sm_75") # Turing
    # set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -gencode arch=compute_86,code=sm_86") # Ampere
    # set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -gencode arch=compute_90,code=sm_90") # Hopper
    # Or target a common minimum:
    # set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -arch=sm_70")

    # --- End TODO CUDA ---

    # Add other desired nvcc flags
    # set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} --use_fast_math") # Be careful with precision

    # Set optimisation level for CUDA based on CMAKE_BUILD_TYPE
    set(CMAKE_CUDA_FLAGS_DEBUG "${CMAKE_CUDA_FLAGS_DEBUG} -G -O0") # Add -G for device debug symbols, -O0 if needed
    set(CMAKE_CUDA_FLAGS_RELEASE "${CMAKE_CUDA_FLAGS_RELEASE} -O3") # Standard release optimization

    message("BUILD CUDA FFLAGS : " )
    message ( ${CMAKE_CUDA_FLAGS} ) # Note: This prints the base flags, not config-specific ones
endif()
# <<< CUDA Change End >>>



### this is for compile_commands.json   (it should help various IDEs and other tool to understand the project)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# ==============================
# ==== Functions
# ==============================

# Attempt to find ccache
find_program(CCACHE_FOUND ccache)
if(CCACHE_FOUND)
    message(STATUS "Found ccache: ${CCACHE_FOUND}")
    # Set ccache as the compiler launcher if found
    set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE "${CCACHE_FOUND}")
else()
    message(STATUS "ccache not found, compilation will proceed without it")
endif()

# ==============================
# ==== DIRS
# ==============================

SET( MY_SRC_DIR    "${CMAKE_SOURCE_DIR}" )
SET( MY_BIN_SHARE  "${MY_SRC_DIR}/Build"       )

SET( COMMON_SRCS      "${MY_SRC_DIR}/common"     )
SET( COMMON_SDL_SRCS  "${MY_SRC_DIR}/common_SDL" )
SET( COMMON_RESOURCES "${MY_SRC_DIR}/common_resources"       )
SET( COMMON_SHADERS   "${MY_SRC_DIR}/common_resources/shaders" )

function(setup_precompiled_headers target)
    set(additional_headers ${ARGN}) # ARGN holds all additional arguments passed to the function
    set(all_headers
        <vector>
        <string>
        <unordered_map>
        <unordered_set>
        "${COMMON_SRCS}/testUtils.h"
        "${COMMON_SRCS}/macroUtils.h"
        "${COMMON_SRCS}/IO_utils.h"
        "${COMMON_SRCS}/math/fastmath.h"
        "${COMMON_SRCS}/math/Vec2.h"
        "${COMMON_SRCS}/math/Vec3.h"
        "${COMMON_SRCS}/math/quaternion.h"
        "${COMMON_SRCS}/math/Mat3.h"
        ${additional_headers} # This now directly uses the list of additional headers
    )
    target_precompile_headers(${target} PRIVATE ${all_headers})
endfunction()

include_directories(
        ${COMMON_SRCS}
        ${COMMON_SRCS}/math
        ${COMMON_SRCS}/OpenCL
        ${COMMON_SRCS}/molecular
        ${COMMON_SRCS}/dataStructures
        ${COMMON_SDL_SRCS}
        ${COMMON_SDL_SRCS}/SDL2
        ${COMMON_SDL_SRCS}/SDL2OGL
)

set( CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/cmake_utils )

#file (GLOB_RECURSE _ALL_SRCS_ *.h *.cpp *.c )
file (GLOB_RECURSE _COMMON_SRCS_     ${COMMON_SRCS}/*.h     ${COMMON_SRCS}/*.cpp     ${COMMON_SRCS}/*.c     )
#file (GLOB_RECURSE _COMMON_SHADERS_  ${COMMON_SHADERS}/*.glslv  ${COMMON_SHADERS}/*.glslf    )

# ==============================
# ==== BUILD TARGETS
#===============================

add_subdirectory( ${MY_SRC_DIR}/common   ${MY_BIN_SHARE}/common    )
add_subdirectory( ${MY_SRC_DIR}/libs    )

# ====== OpenCL

if ( WITH_OPENCL )
    add_subdirectory( ${MY_SRC_DIR}/libs_OCL )
    find_package( OpenCL REQUIRED )
endif()

if ( WITH_CUDA )
    add_subdirectory( ${MY_SRC_DIR}/apps_CUDA )
endif()





# ====== SDL@OpenGL

if ( WITH_SDL )
    # https://cmake.org/cmake/help/latest/policy/CMP0072.html
    #set(OpenGL_GL_PREFERENCE LEGACY)
    set(OpenGL_GL_PREFERENCE GLVND)
    file (GLOB_RECURSE _COMMON_SDL_SRCS_ ${COMMON_SDL_SRCS}/*.h ${COMMON_SDL_SRCS}/*.cpp ${COMMON_SDL_SRCS}/*.c )
    add_subdirectory( ${MY_SRC_DIR}/common_SDL   ${MY_BIN_SHARE}/common_SDL )
    add_subdirectory( ${MY_SRC_DIR}/libs_SDL     )
    add_subdirectory( ${MY_SRC_DIR}/sketches_SDL )
    add_subdirectory( ${MY_SRC_DIR}/apps         )
endif()

# ====== OpenCL + SDL@OpenGL

if ( WITH_OPENCL AND WITH_SDL )
    add_subdirectory( ${MY_SRC_DIR}/apps_OCL     )
endif()