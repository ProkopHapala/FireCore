
project(Apps)

find_package( OpenGL REQUIRED )
find_package( SDL2   REQUIRED )

add_subdirectory (MolecularEditor)
add_subdirectory (EFF)
add_subdirectory (RFF)



file (GLOB_RECURSE _COMMON_SRCS_     ${COMMON_SRCS}/*.h     ${COMMON_SRCS}/*.cpp     ${COMMON_SRCS}/*.c     )
file (GLOB_RECURSE _COMMON_SDL_SRCS_ ${COMMON_SDL_SRCS}/*.h ${COMMON_SDL_SRCS}/*.cpp ${COMMON_SDL_SRCS}/*.c )

add_custom_target( ${PROJECT_NAME}_src_dep SOURCES ${_COMMON_SRCS_} ${_COMMON_SDL_SRCS_} )
