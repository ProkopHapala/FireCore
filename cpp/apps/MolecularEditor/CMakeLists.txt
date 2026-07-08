
set( NAME "MolecularEditor" )

#include_directories( ${COMMON_SRCS}/dynamics/Molecular )
include_directories( ${COMMON_SRCS}/molecular )
include_directories( ${COMMON_SRCS}/io )

add_executable( MolecularBrowser
    MolecularBrowser.cpp 
	$<TARGET_OBJECTS:DynamicOpt>
	$<TARGET_OBJECTS:SDL2OGL>  
)
setup_precompiled_headers( MolecularBrowser )
target_link_libraries( MolecularBrowser ${OPENGL_LIBRARIES} ${SDL2_LIBRARIES} ${LUA_LIBRARIES} z )

add_executable( MolGUIapp
	MolGUIapp.cpp
	$<TARGET_OBJECTS:DynamicOpt>
	$<TARGET_OBJECTS:SDL2OGL>
    $<TARGET_OBJECTS:ProjectiveDynamics_d>    
)
setup_precompiled_headers( MolGUIapp )
target_link_libraries( MolGUIapp ${OPENGL_LIBRARIES} ${SDL2_LIBRARIES} ${LUA_LIBRARIES} ${FFTW_LIBRARIES} )


add_executable( MolGUIapp_QMMM
	MolGUIapp_QMMM.cpp
	$<TARGET_OBJECTS:DynamicOpt>
	$<TARGET_OBJECTS:SDL2OGL>  
    $<TARGET_OBJECTS:ProjectiveDynamics_d> 
)
setup_precompiled_headers( MolGUIapp_QMMM )
target_link_libraries( MolGUIapp_QMMM ${OPENGL_LIBRARIES} ${SDL2_LIBRARIES} ${LUA_LIBRARIES} ${FFTW_LIBRARIES} )

#find_package(GLEW REQUIRED)
#if (GLEW_FOUND)
#    include_directories(${GLEW_INCLUDE_DIRS})
#    link_libraries(${GLEW_LIBRARIES})
#endif()
#add_executable       ( test_CLCFGO_OGL3 test_CLCFGO_OGL3.cpp $<TARGET_OBJECTS:SDL2OGL3> )
#target_link_libraries( test_CLCFGO_OGL3 ${OPENGL_LIBRARIES}   ${GLEW_LIBRARY} ${SDL2_LIBRARIES} )
#target_link_libraries( test_CLCFGO_OGL3 ${OPENGL_LIBRARIES} ${GLEW_LIBRARY} ${SDL2_LIBRARIES} )

add_custom_command(
    TARGET ${NAME}_main ${NAME}_old_main ConfSearch test_MMFF MolecularBrowser
    COMMAND ${CMAKE_COMMAND} -E create_symlink ${MY_SRC_DIR}/common_resources     ${CMAKE_CURRENT_BINARY_DIR}/common_resources
    COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_CURRENT_SOURCE_DIR}/data   ${CMAKE_CURRENT_BINARY_DIR}/data
    COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_CURRENT_SOURCE_DIR}/inputs ${CMAKE_CURRENT_BINARY_DIR}/inputs
)
