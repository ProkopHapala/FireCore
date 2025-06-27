set( NAME "apps_CUDA" )

find_package(CUDA REQUIRED)
include_directories(${CUDA_INCLUDE_DIRS})
include_directories(${COMMON_RESOURCES}/cl)

#add_library( MolecularEngine OBJECT MolecularWorld.cpp MoleculeType.cpp AtomTypes.cpp  )
SET( MOLECULAR_EDITOR_SRCS       ${MY_SRC_DIR}/apps/MolecularEditor )

include_directories( ${COMMON_SRCS}/molecular ${MOLECULAR_EDITOR_SRCS} )
#include_directories( ${COMMON_SRCS}/dynamics/Molecular ${MOLECULAR_EDITOR_SRCS} )

# Rename executable if desired
# add_executable( MolGUIapp_cuda # Renamed
#         MolGUIapp_cuda.cpp # Keep if common, otherwise rename/replace
#         # --- TODO CUDA: List all C++, CUDA (.cu, .cuh) source files ---
#         $<TARGET_OBJECTS:DynamicOpt>
#         $<TARGET_OBJECTS:SDL2OGL>
#         $<TARGET_OBJECTS:ProjectiveDynamics_d>
# )
# # Link libraries: Include CUDA runtime
# target_link_libraries( MolGUIapp_cuda # Renamed target
#     ${CUDA_LIBRARIES}
#     # PUBLIC CUDA::cudart # Example explicit link if needed
#     # PUBLIC CUDA::cufft # Example explicit link if needed
#     ${OPENGL_LIBRARIES}
#     # ${CLFFT_LIBRARIES} # Remove or replace
#     ${SDL2_LIBRARIES}
#     ${LUA_LIBRARIES}
#     ${FFTW_LIBRARIES}
# )

add_library(cuMMFF_lib  SHARED cuMMFF_lib.cpp ${COMMON_RESOURCES}/cl/relax_multi.cu)
set_source_files_properties(
    cuMMFF_lib.cpp
    PROPERTIES LANGUAGE CUDA
)
set_source_files_properties(${COMMON_RESOURCES}/cl/relax_multi.cu PROPERTIES LANGUAGE CUDA)
set_target_properties(cuMMFF_lib PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
target_link_libraries(cuMMFF_lib ${CUDA_LIBRARIES})

add_custom_command(
    TARGET cuMMFF_lib
#    TARGET MolGUIapp_cuda
#    TARGET ${NAME}_main ${NAME}_scanner
    COMMAND ${CMAKE_COMMAND} -E create_symlink ${MY_SRC_DIR}/common_resources  ${CMAKE_CURRENT_BINARY_DIR}/common_resources
    COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_CURRENT_SOURCE_DIR}/cl  ${CMAKE_CURRENT_BINARY_DIR}/cl
    COMMAND ${CMAKE_COMMAND} -E create_symlink ${MOLECULAR_EDITOR_SRCS}/inputs ${CMAKE_CURRENT_BINARY_DIR}/inputs
    COMMAND ${CMAKE_COMMAND} -E create_symlink ${MOLECULAR_EDITOR_SRCS}/inputs ${CMAKE_CURRENT_SOURCE_DIR}/inputs
    COMMAND ${CMAKE_COMMAND} -E create_symlink ${MOLECULAR_EDITOR_SRCS}/data   ${CMAKE_CURRENT_BINARY_DIR}/data
    COMMAND ${CMAKE_COMMAND} -E create_symlink ${MOLECULAR_EDITOR_SRCS}/data   ${CMAKE_CURRENT_SOURCE_DIR}/data
    COMMAND ${CMAKE_COMMAND} -E create_symlink ${MOLECULAR_EDITOR_SRCS}/doc    ${CMAKE_CURRENT_BINARY_DIR}/doc
    COMMAND ${CMAKE_COMMAND} -E create_symlink ${MOLECULAR_EDITOR_SRCS}/doc    ${CMAKE_CURRENT_SOURCE_DIR}/doc
)
