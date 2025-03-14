
set( NAME "EFF" )

#include_directories( ${COMMON_SRCS}/dynamics/Molecular )
include_directories( ${COMMON_SRCS}/molecular )

add_executable( EFFapp EFFapp.cpp $<TARGET_OBJECTS:DynamicOpt> $<TARGET_OBJECTS:SDL2OGL> )
target_link_libraries( EFFapp ${OPENGL_LIBRARIES} ${SDL2_LIBRARIES} ${LUA_LIBRARIES} )

#add_executable( EFFapp_console EFFapp_console.cpp $<TARGET_OBJECTS:DynamicOpt> )
#target_link_libraries( EFFapp_console ${LUA_LIBRARIES} )


add_library(EFFapp_console SHARED EFFapp_console.cpp $<TARGET_OBJECTS:DynamicOpt>)
target_link_libraries(EFFapp_console ${LUA_LIBRARIES})
set_target_properties(EFFapp_console PROPERTIES 
    LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/apps/EFF
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/apps/EFF
)


add_custom_command(
    TARGET EFFapp
    COMMAND ${CMAKE_COMMAND} -E create_symlink ${MY_SRC_DIR}/common_resources             ${CMAKE_CURRENT_BINARY_DIR}/common_resources
    COMMAND ${CMAKE_COMMAND} -E create_symlink ${MY_SRC_DIR}/sketches_SDL/Molecular/data  ${CMAKE_CURRENT_BINARY_DIR}/data
    COMMAND ${CMAKE_COMMAND} -E create_symlink ${MY_SRC_DIR}/sketches_SDL/Molecular/data  ${CMAKE_CURRENT_BINARY_DIR}/inputs
)
