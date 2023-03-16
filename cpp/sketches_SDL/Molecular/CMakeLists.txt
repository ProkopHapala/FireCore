
project(sketches_SDL_Molecular)

#include_directories( ${COMMON_SRCS}/dynamics/Molecular )
include_directories( ${COMMON_SRCS}/molecular )

add_executable       ( test_sp3space test_sp3space.cpp $<TARGET_OBJECTS:SDL2OGL> $<TARGET_OBJECTS:DynamicOpt>  )
target_link_libraries( test_sp3space ${OPENGL_LIBRARIES}   ${SDL2_LIBRARIES} )

add_executable       ( test_spRotations test_spRotations.cpp $<TARGET_OBJECTS:SDL2OGL> $<TARGET_OBJECTS:DynamicOpt> $<TARGET_OBJECTS:Lingebra>  )
target_link_libraries( test_spRotations ${OPENGL_LIBRARIES}   ${SDL2_LIBRARIES} )

add_executable       ( test_CLCFSF test_CLCFSF.cpp $<TARGET_OBJECTS:SDL2OGL> $<TARGET_OBJECTS:DynamicOpt> $<TARGET_OBJECTS:Lingebra>  )
target_link_libraries( test_CLCFSF ${OPENGL_LIBRARIES}   ${SDL2_LIBRARIES} )

add_executable       ( test_CLCFGO test_CLCFGO.cpp $<TARGET_OBJECTS:SDL2OGL> $<TARGET_OBJECTS:DynamicOpt> $<TARGET_OBJECTS:Lingebra>  )
target_link_libraries( test_CLCFGO ${OPENGL_LIBRARIES}   ${SDL2_LIBRARIES} )

add_executable       ( test_ConfDynamics test_ConfDynamics.cpp $<TARGET_OBJECTS:SDL2OGL> $<TARGET_OBJECTS:DynamicOpt> )
target_link_libraries( test_ConfDynamics ${OPENGL_LIBRARIES}   ${SDL2_LIBRARIES} )

add_executable       ( test_Multipoles test_Multipoles.cpp $<TARGET_OBJECTS:SDL2OGL> )
target_link_libraries( test_Multipoles ${OPENGL_LIBRARIES}   ${SDL2_LIBRARIES} )

add_executable       ( test_SoftMolecularDynamics test_SoftMolecularDynamics.cpp $<TARGET_OBJECTS:SDL2OGL> $<TARGET_OBJECTS:DynamicOpt>  )
target_link_libraries( test_SoftMolecularDynamics ${OPENGL_LIBRARIES}   ${SDL2_LIBRARIES} )

add_executable       ( test_FIRE test_FIRE.cpp $<TARGET_OBJECTS:SDL2OGL> $<TARGET_OBJECTS:DynamicOpt>  )
target_link_libraries( test_FIRE ${OPENGL_LIBRARIES}   ${SDL2_LIBRARIES} )

add_executable       ( test_eFF test_eFF.cpp $<TARGET_OBJECTS:SDL2OGL> $<TARGET_OBJECTS:DynamicOpt>  )
target_link_libraries( test_eFF ${OPENGL_LIBRARIES}   ${SDL2_LIBRARIES} )

add_executable       ( test_eFF_old test_eFF_old.cpp $<TARGET_OBJECTS:SDL2OGL> $<TARGET_OBJECTS:DynamicOpt>  )
target_link_libraries( test_eFF_old ${OPENGL_LIBRARIES}   ${SDL2_LIBRARIES} )

add_executable       ( test_eFF_MD test_eFF_MD.cpp $<TARGET_OBJECTS:SDL2OGL> $<TARGET_OBJECTS:DynamicOpt>  )
target_link_libraries( test_eFF_MD ${OPENGL_LIBRARIES}   ${SDL2_LIBRARIES} )

add_executable       ( test_eFFMC test_eFFMC.cpp $<TARGET_OBJECTS:SDL2OGL> $<TARGET_OBJECTS:DynamicOpt>  )
target_link_libraries( test_eFFMC ${OPENGL_LIBRARIES}   ${SDL2_LIBRARIES} )

add_executable       ( test_MMFFmini test_MMFFmini.cpp $<TARGET_OBJECTS:SDL2OGL> $<TARGET_OBJECTS:DynamicOpt>  )
target_link_libraries( test_MMFFmini ${OPENGL_LIBRARIES}   ${SDL2_LIBRARIES} )

add_executable       ( test_MMFFsp3 test_MMFFsp3.cpp $<TARGET_OBJECTS:SDL2OGL> $<TARGET_OBJECTS:DynamicOpt>  )
target_link_libraries( test_MMFFsp3 ${OPENGL_LIBRARIES}   ${SDL2_LIBRARIES} )

add_executable       ( test_FF2D test_FF2D.cpp $<TARGET_OBJECTS:SDL2OGL> $<TARGET_OBJECTS:DynamicOpt>  )
target_link_libraries( test_FF2D ${OPENGL_LIBRARIES}   ${SDL2_LIBRARIES} )

add_executable       ( test_DirectionStiffness test_DirectionStiffness.cpp $<TARGET_OBJECTS:SDL2OGL> $<TARGET_OBJECTS:DynamicOpt>  )
target_link_libraries(test_DirectionStiffness  ${OPENGL_LIBRARIES}   ${SDL2_LIBRARIES} )

add_executable       ( test_FARFF test_FARFF.cpp $<TARGET_OBJECTS:SDL2OGL> $<TARGET_OBJECTS:DynamicOpt>  )
target_link_libraries( test_FARFF ${OPENGL_LIBRARIES}   ${SDL2_LIBRARIES} )

add_executable       ( test_RARFF test_RARFF.cpp $<TARGET_OBJECTS:SDL2OGL> $<TARGET_OBJECTS:DynamicOpt>  )
target_link_libraries( test_RARFF ${OPENGL_LIBRARIES}   ${SDL2_LIBRARIES} )

add_executable       ( test_RARFF_SR test_RARFF_SR.cpp $<TARGET_OBJECTS:SDL2OGL> $<TARGET_OBJECTS:DynamicOpt>  )
target_link_libraries( test_RARFF_SR ${OPENGL_LIBRARIES}   ${SDL2_LIBRARIES} )

add_executable       ( test_RARFFarr test_RARFFarr.cpp $<TARGET_OBJECTS:SDL2OGL> $<TARGET_OBJECTS:DynamicOpt>  )
target_link_libraries( test_RARFFarr ${OPENGL_LIBRARIES}   ${SDL2_LIBRARIES} )

add_executable       ( test_RARFF2 test_RARFF2.cpp $<TARGET_OBJECTS:SDL2OGL> $<TARGET_OBJECTS:DynamicOpt>  )
target_link_libraries( test_RARFF2 ${OPENGL_LIBRARIES}   ${SDL2_LIBRARIES} )

add_executable       ( test_RspFF test_RspFF.cpp $<TARGET_OBJECTS:SDL2OGL> $<TARGET_OBJECTS:DynamicOpt>  )
target_link_libraries( test_RspFF ${OPENGL_LIBRARIES}   ${SDL2_LIBRARIES} )

add_executable       ( test_FTRFF test_FTRFF.cpp $<TARGET_OBJECTS:SDL2OGL> $<TARGET_OBJECTS:DynamicOpt>  )
target_link_libraries( test_FTRFF ${OPENGL_LIBRARIES}   ${SDL2_LIBRARIES} )

add_executable       ( test_RRFF test_RRFF.cpp $<TARGET_OBJECTS:SDL2OGL> $<TARGET_OBJECTS:DynamicOpt>  )
target_link_libraries( test_RRFF ${OPENGL_LIBRARIES}   ${SDL2_LIBRARIES} )

add_executable       ( test_EOFF test_EOFF.cpp $<TARGET_OBJECTS:SDL2OGL> $<TARGET_OBJECTS:DynamicOpt>  )
target_link_libraries( test_EOFF ${OPENGL_LIBRARIES}   ${SDL2_LIBRARIES} )

add_executable       ( test_ESFF test_ESFF.cpp $<TARGET_OBJECTS:SDL2OGL> $<TARGET_OBJECTS:DynamicOpt>  )
target_link_libraries( test_ESFF ${OPENGL_LIBRARIES}   ${SDL2_LIBRARIES} )

add_executable       ( test_BondAdaptedMesh test_BondAdaptedMesh.cpp $<TARGET_OBJECTS:SDL2OGL>  )
target_link_libraries( test_BondAdaptedMesh ${OPENGL_LIBRARIES}   ${SDL2_LIBRARIES} )

add_executable       ( test_MatchLattice2D test_MatchLattice2D.cpp $<TARGET_OBJECTS:SDL2OGL> )
target_link_libraries( test_MatchLattice2D ${OPENGL_LIBRARIES}   ${SDL2_LIBRARIES} )

add_custom_target( ${PROJECT_NAME}_src_dep SOURCES ${_COMMON_SRCS_} ${_COMMON_SDL_SRCS_} )

add_custom_command(
    TARGET ${PROJECT_NAME} test_SoftMolecularDynamics test_eFF test_eFF_old
    #COMMAND ln -sf ${CMAKE_SOURCE_DIR}/common_resources ${CMAKE_CURRENT_BINARY_DIR}/common_resources
    #COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_SOURCE_DIR}/common_resources ${CMAKE_CURRENT_BINARY_DIR}/common_resources
    COMMAND ${CMAKE_COMMAND} -E create_symlink ${MY_SRC_DIR}/common_resources ${CMAKE_CURRENT_BINARY_DIR}/common_resources
    COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_CURRENT_SOURCE_DIR}/data ${CMAKE_CURRENT_BINARY_DIR}/data
    COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_CURRENT_BINARY_DIR}/temp
)
