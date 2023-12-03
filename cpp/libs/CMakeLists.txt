
# problem is that these libraries cannot be build from common objects .o; needs "-fPIC" to relocation

project(libs)

#add_subdirectory (libSail)
#add_subdirectory (libFlight)
add_subdirectory (Molecular)
#add_subdirectory (Shock)
#add_subdirectory (KosmoSuite)
#add_subdirectory (CombatModels)

#add_library(quadrature_lib  SHARED quadrature_lib.cpp  $<TARGET_OBJECTS:DynamicOpt> )
add_library(quadrature_lib  SHARED quadrature_lib.cpp  )
