
set( NAME "MolecularEditorOCL" )

#add_library( MolecularEngine OBJECT MolecularWorld.cpp MoleculeType.cpp AtomTypes.cpp  )
SET( MOLECULAR_EDITOR_SRCS       ${MY_SRC_DIR}/apps/MolecularEditor )

include_directories( ${COMMON_SRCS}/molecular ${MOLECULAR_EDITOR_SRCS} )
#include_directories( ${COMMON_SRCS}/dynamics/Molecular ${MOLECULAR_EDITOR_SRCS} )

add_executable( FireCoreVisualOCL
FireCoreVisualOCL.cpp
	$<TARGET_OBJECTS:DynamicOpt>
	$<TARGET_OBJECTS:SDL2OGL>  
)
target_link_libraries( FireCoreVisualOCL  ${OPENCL_LIBRARY}  ${CLFFT_LIBRARIES} ${OPENGL_LIBRARIES}   ${SDL2_LIBRARIES} ${LUA_LIBRARIES} -rdynamic -ldl -lm )

add_executable( MolGUIapp_multi
        MolGUIapp_multi.cpp
        $<TARGET_OBJECTS:DynamicOpt>
        $<TARGET_OBJECTS:SDL2OGL>  
)
target_link_libraries( MolGUIapp_multi ${OPENGL_LIBRARIES} ${CLFFT_LIBRARIES} ${OPENGL_LIBRARIES} ${SDL2_LIBRARIES} ${LUA_LIBRARIES} ${FFTW_LIBRARIES} )

add_executable( MolGUIapp_QMMM_multi
        MolGUIapp_QMMM_multi.cpp
        $<TARGET_OBJECTS:DynamicOpt>
        $<TARGET_OBJECTS:SDL2OGL>  
)
target_link_libraries( MolGUIapp_QMMM_multi ${OPENGL_LIBRARIES} ${CLFFT_LIBRARIES} ${OPENGL_LIBRARIES} ${SDL2_LIBRARIES} ${LUA_LIBRARIES} ${FFTW_LIBRARIES} )

#add_custom_target( ${NAME}_src_dep SOURCES ${MOLECULAR_EDITOR_SRCS} )

add_custom_command(
    TARGET FireCoreVisualOCL MolGUIapp_multi MolGUIapp_QMMM_multi
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
