
set( NAME "RFF" )

#include_directories( ${COMMON_SRCS}/dynamics/Molecular )
include_directories( ${COMMON_SRCS}/molecular )

add_executable( RFFapp.x
    RFFapp.cpp 
	$<TARGET_OBJECTS:DynamicOpt>
	$<TARGET_OBJECTS:SDL2OGL>  
)
target_link_libraries( RFFapp.x ${OPENGL_LIBRARIES} ${SDL2_LIBRARIES} ${LUA_LIBRARIES} )

#add_executable       ( test_eFF test_eFF.cpp $<TARGET_OBJECTS:SDL2OGL> $<TARGET_OBJECTS:DynamicOpt>  )
#target_link_libraries( test_eFF ${OpenGL_LIBRARY} ${GLU_LIBRARY} ${SDL2_LIBRARY} )

#add_executable       ( test_RARFF_SR test_RARFF_SR.cpp $<TARGET_OBJECTS:SDL2OGL> $<TARGET_OBJECTS:DynamicOpt>  )
#target_link_libraries( test_RARFF_SR ${OpenGL_LIBRARY} ${GLU_LIBRARY} ${SDL2_LIBRARY} )

add_custom_command(
    TARGET RFFapp.x
    COMMAND ${CMAKE_COMMAND} -E create_symlink ${MY_SRC_DIR}/common_resources             ${CMAKE_CURRENT_BINARY_DIR}/common_resources
    COMMAND ${CMAKE_COMMAND} -E create_symlink ${MY_SRC_DIR}/sketches_SDL/Molecular/data  ${CMAKE_CURRENT_BINARY_DIR}/data
    COMMAND ${CMAKE_COMMAND} -E create_symlink ${MY_SRC_DIR}/sketches_SDL/Molecular/data  ${CMAKE_CURRENT_BINARY_DIR}/inputs
)
