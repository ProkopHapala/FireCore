
set( NAME "MolecularEditor" )

#include_directories( ${COMMON_SRCS}/dynamics/Molecular )
include_directories( ${COMMON_SRCS}/molecular )

add_executable( ${NAME}_main
	${NAME}_main.cpp 
	$<TARGET_OBJECTS:DynamicOpt>
	$<TARGET_OBJECTS:SDL2OGL>  
)
setup_precompiled_headers(${NAME}_main )
target_link_libraries( ${NAME}_main ${OPENGL_LIBRARIES} ${SDL2_LIBRARIES} ${LUA_LIBRARIES} )


add_executable( MolecularBrowser
    MolecularBrowser.cpp 
	$<TARGET_OBJECTS:DynamicOpt>
	$<TARGET_OBJECTS:SDL2OGL>  
)
setup_precompiled_headers( MolecularBrowser )
target_link_libraries( MolecularBrowser ${OPENGL_LIBRARIES} ${SDL2_LIBRARIES} ${LUA_LIBRARIES} )


add_executable( FireCoreVisual
	FireCoreVisual.cpp
	$<TARGET_OBJECTS:DynamicOpt>
	$<TARGET_OBJECTS:SDL2OGL>  
)
setup_precompiled_headers( FireCoreVisual )
target_link_libraries( FireCoreVisual ${OPENGL_LIBRARIES} ${SDL2_LIBRARIES} ${LUA_LIBRARIES} -rdynamic -ldl -lm )


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


add_executable( FireCoreMMFFmini
	FireCoreMMFFmini.cpp
	$<TARGET_OBJECTS:DynamicOpt>
	$<TARGET_OBJECTS:SDL2OGL>  
)
setup_precompiled_headers( FireCoreMMFFmini )
target_link_libraries( FireCoreMMFFmini ${OPENGL_LIBRARIES} ${SDL2_LIBRARIES} ${LUA_LIBRARIES} -rdynamic -ldl -lm )


#add_executable( ${NAME}_old_main
#	${NAME}_old_main.cpp 
#	#$<TARGET_OBJECTS:Body> 
#	$<TARGET_OBJECTS:MolecularEngine>  
#	$<TARGET_OBJECTS:DynamicOpt>
#	$<TARGET_OBJECTS:SDL2OGL>  
#)
#setup_precompiled_headers(  )
#target_link_libraries( ${NAME}_old_main ${OPENGL_LIBRARIES} ${SDL2_LIBRARIES} )


add_executable( ConfSearch 
	ConfSearch.cpp 
	$<TARGET_OBJECTS:DynamicOpt>
	$<TARGET_OBJECTS:SDL2OGL>  
)
setup_precompiled_headers( ConfSearch )
target_link_libraries( ConfSearch ${OPENGL_LIBRARIES} ${SDL2_LIBRARIES} ${LUA_LIBRARIES} )

#add_executable( test_reset 	test_reset.cpp 	$<TARGET_OBJECTS:DynamicOpt> 	$<TARGET_OBJECTS:SDL2OGL>   )
#target_link_libraries( test_reset ${OPENGL_LIBRARIES}   ${SDL2_LIBRARIES} )

add_executable( test_MMFF  test_MMFF.cpp  $<TARGET_OBJECTS:DynamicOpt> )
setup_precompiled_headers( test_MMFF )
target_link_libraries( test_MMFF )

#find_package(GLEW REQUIRED)
#if (GLEW_FOUND)
#    include_directories(${GLEW_INCLUDE_DIRS})
#    link_libraries(${GLEW_LIBRARIES})
#endif()
#add_executable       ( test_CLCFGO_OGL3 test_CLCFGO_OGL3.cpp $<TARGET_OBJECTS:SDL2OGL3> )
#target_link_libraries( test_CLCFGO_OGL3 ${OPENGL_LIBRARIES}   ${GLEW_LIBRARY} ${SDL2_LIBRARIES} )
#target_link_libraries( test_CLCFGO_OGL3 ${OPENGL_LIBRARIES} ${GLEW_LIBRARY} ${SDL2_LIBRARIES} )

add_custom_command(
    TARGET ${NAME}_main ${NAME}_old_main ConfSearch test_MMFF
    COMMAND ${CMAKE_COMMAND} -E create_symlink ${MY_SRC_DIR}/common_resources     ${CMAKE_CURRENT_BINARY_DIR}/common_resources
    COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_CURRENT_SOURCE_DIR}/data   ${CMAKE_CURRENT_BINARY_DIR}/data
    COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_CURRENT_SOURCE_DIR}/inputs ${CMAKE_CURRENT_BINARY_DIR}/inputs
)
