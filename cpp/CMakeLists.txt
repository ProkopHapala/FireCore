#cmake_minimum_required ( VERSION 2.8 )
cmake_minimum_required(VERSION 3.10)

project ( SimpleSimulationEngine )

option(WITH_SDL    "build with GUI and 3D graphics (SDL+OpenGL)"                           OFF)
option(WITH_LUA    "build with Lua 5.2 for scripting inside interactive GUI"               OFF)
option(WITH_FFTW   "build with FFTW3 for scripting inside interactive GUI"                 OFF)
option(WITH_OPENCL "build with OpenCL GPU acceleration? up to 100x speedup"                OFF)
option(WITH_OMP    "build with OpenMP paralelization for multiple CPUs"                    OFF)
option(WITH_AVX    "build with AVX256 SIMD intrinstic ? up to 10x speedup on 1CPU"         OFF)
option(RELEASE     "build optimized release (-O3/-Ofast) rather than debug? (-g -Og)"      OFF)
option(WITH_ASAN   "use runtime memory sanitizer (asan) to trace segfaults and mem-leaks?" OFF)
option(TARGET_WASM "compile to WASM for running in html5."                                 OFF)

message("OPTIONS: -DWITH_SDL=${WITH_SDL} -DWITH_LUA=${WITH_LUA} -DWITH_FFTW=${WITH_FFTW} -DWITH_OPENCL=${WITH_OPENCL} -DWITH_OMP=${WITH_OMP} -DWITH_AVX=${WITH_AVX} -DRELEASE=${RELEASE} -DWITH_ASAN=${WITH_ASAN} -DTARGET_WASM=${TARGET_WASM}" )

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

if(TARGET_WASM)
	if(WITH_FFTW)
		message( FATAL_ERROR "-DWITH_FFTW must be 'NO' when targeting Web Assembly (-DTARGET_WASM=ON)." )
	endif()
	if(WITH_OMP)
		message( FATAL_ERROR "-DWITH_OMP must be 'NO' when targeting Web Assembly (-DTARGET_WASM=ON)." )
	endif()

	SET(CMAKE_EXECUTABLE_SUFFIX ".html")
	SET(AXULIARY_COMPILE_FLAGS "${AXULIARY_COMPILE_FLAGS} -fno-openmp -std=c++20 -march=wasm")
	SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} --use-port=sdl2 -sUSE_SDL=2 -sMIN_WEBGL_VERSION=3 -sMAX_WEBGL_VERSION=3 -sALLOW_MEMORY_GROWTH=1")
	SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} --embed-file ${CMAKE_SOURCE_DIR}/common_resources@common_resources --embed-file ${CMAKE_SOURCE_DIR}/common_resources@data")
	SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -sTOTAL_MEMORY=268435456") # 256MiB
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
	SET( AXULIARY_COMPILE_FLAGS "${AXULIARY_COMPILE_FLAGS} -g")
	if ( TARGET_WASM ) # emscripten doesn't work with -Og, and needs -march=wasm
		SET( AXULIARY_COMPILE_FLAGS "${AXULIARY_COMPILE_FLAGS} -O0")
	else()
		SET( AXULIARY_COMPILE_FLAGS "${AXULIARY_COMPILE_FLAGS} -Og -march=native")
	endif()
	
    if ( WITH_ASAN )
        SET( AXULIARY_COMPILE_FLAGS "${AXULIARY_COMPILE_FLAGS} -fsanitize=address -fno-omit-frame-pointer")
        #SET( AXULIARY_COMPILE_FLAGS "-g -Og -fsanitize=address -fno-omit-frame-pointer -static-libasan") # static-asan does not work see https://github.com/google/sanitizers/issues/796
    endif()
    #set(CMAKE_CXX_CPPCHECK "cppcheck")
    #message("DEBUG BUILD FFLAGS : " )
else()
    #SET( AXULIARY_COMPILE_FLAGS "-w -O3 -fPIC ")
    #SET( AXULIARY_COMPILE_FLAGS "-w -O2 -march=native")
    #SET( AXULIARY_COMPILE_FLAGS "-w -O2 -march=native")
    SET( AXULIARY_COMPILE_FLAGS "${AXULIARY_COMPILE_FLAGS} -w -Ofast -march=native -mtune=native")
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
SET( IGNORE_WARRNING_FLAGS "-Wno-int-to-pointer-cast -Wno-switch -Wno-reorder -Wno-deprecated -Wno-char-subscripts -Wno-write-strings -Wno-format -Wno-parentheses -Wno-unused-but-set-variable  -Wno-narrowing -Wno-unused-result  -Wno-sign-compare -Wno-strict-aliasing -Wno-unused-variable -Wno-unused-value -Wno-comment -Wno-misleading-indentation " )
if ( TARGET_WASM )
	SET( IGNORE_WARRNING_FLAGS "${IGNORE_WARRNING_FLAGS} -Wno-implicit-const-int-float-conversion" )
else ()
	SET( IGNORE_WARRNING_FLAGS "${IGNORE_WARRNING_FLAGS} -Wno-maybe-uninitialized" )
endif()
SET( WARRNING_TO_ERROR     "-Werror=return-type -Werror=init-self -Werror=uninitialized -Werror=implicit-fallthrough -Werror=tautological-compare -Werror=delete-non-virtual-dtor")
SET( CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${AXULIARY_COMPILE_FLAGS} ${IGNORE_WARRNING_FLAGS} ${WARRNING_TO_ERROR}" )
message("BUILD FFLAGS : " )
message ( ${CMAKE_CXX_FLAGS} )

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

if ( TARGET_WASM )
include_directories(
        ${COMMON_SRCS}
        ${COMMON_SRCS}/math
        ${COMMON_SRCS}/OpenCL
        ${COMMON_SRCS}/molecular
        ${COMMON_SRCS}/dataStructures
        ${COMMON_SDL_SRCS}
        ${COMMON_SDL_SRCS}/SDL2
        ${COMMON_SDL_SRCS}/SDL2OGL
        
        ${CMAKE_CURRENT_BINARY_DIR}/simde
)
else()
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
endif()

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

if(TARGET_WASM)
	ADD_CUSTOM_COMMAND(OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/simde/simde/x86/fma.h"
		VERBATIM COMMAND cd ${CMAKE_CURRENT_BINARY_DIR} && git clone https://github.com/simd-everywhere/simde 
	)
	ADD_CUSTOM_TARGET(simde DEPENDS "${CMAKE_CURRENT_BINARY_DIR}/simde/simde/x86/fma.h")
	add_dependencies(MolGUIapp simde)

	#ADD_CUSTOM_TARGET(vfilesystem
	#	VERBATIM COMMAND "/usr/share/emscripten/tools/file_packager TARGET --embed ${CMAKE_SOURCE_DIR}/common_resources@common_resources --embed ${CMAKE_SOURCE_DIR}/common_resources@data --js-output=${CMAKE_CURRENT_BINARY_DIR}/vfilesystem.js"
	#)
	#add_dependencies(MolGUIapp vfilesystem)

	## bugfix: not sure why it tries to link -lnul, but it doesn't exist and we don't need it
	get_target_property(TARGET_LIBRARIES MolGUIapp LINK_LIBRARIES)
	LIST(REMOVE_ITEM TARGET_LIBRARIES nul )
	set_property(TARGET MolGUIapp PROPERTY LINK_LIBRARIES  ${TARGET_LIBRARIES} )
endif()
