
#include_directories( ${COMMON_SRCS}/dynamics/Molecular )
include_directories( ${COMMON_SRCS}/molecular )
#add_library(Molecular   SHARED Molecular.cpp  $<TARGET_OBJECTS:MolecularEngine> $<TARGET_OBJECTS:DynamicOpt> )
#add_library(RigidMol    SHARED RigidMol.cpp   $<TARGET_OBJECTS:MolecularEngine> $<TARGET_OBJECTS:DynamicOpt> )
add_library(ReactiveFF  SHARED ReactiveFF.cpp )

add_library(CLCFGO_lib   SHARED CLCFGO_lib.cpp )
add_library(eFF_lib      SHARED eFF_lib.cpp      $<TARGET_OBJECTS:DynamicOpt> )

# ---- Tomas's original (pre-port) eFF build; pristine sources kept in origZ/ for A/B testing (see PORT_eFF_optimization.md) ----
# origZ dirs MUST precede the global common/ includes so the original fastmath.h / DynamicOpt.h / eFF.h shadow the ported ones
set(ORIGZ ${CMAKE_CURRENT_SOURCE_DIR}/origZ)
add_library( DynamicOpt_origZ OBJECT ${ORIGZ}/math/DynamicOpt.cpp )
target_include_directories( DynamicOpt_origZ BEFORE PRIVATE ${ORIGZ}/math ${ORIGZ}/molecular )
add_library( eFF_lib_origZ SHARED ${ORIGZ}/eFF_lib.cpp $<TARGET_OBJECTS:DynamicOpt_origZ> )
target_include_directories( eFF_lib_origZ BEFORE PRIVATE ${ORIGZ}/math ${ORIGZ}/molecular )
add_library(MMFF_lib     SHARED MMFF_lib.cpp     $<TARGET_OBJECTS:DynamicOpt> $<TARGET_OBJECTS:ProjectiveDynamics_d> )
target_link_libraries(MMFF_lib ${FFTW_LIBRARIES})

#add_library(MMFF_lib     SHARED MMFF_lib.cpp     $<TARGET_OBJECTS:DynamicOpt> )
add_library(MMFFsp3_lib  SHARED MMFFsp3_lib.cpp  $<TARGET_OBJECTS:DynamicOpt> )

add_library(FF2D_lib                 SHARED FF2D_lib.cpp                    )
add_library(Kekule_lib               SHARED Kekule_lib.cpp              )
add_library(Lattice2D_lib            SHARED Lattice2D_lib.cpp           )
add_library(SchroedingerGreen2D_lib  SHARED SchroedingerGreen2D_lib.cpp )
add_library(SchroedingerGreen1D_lib  SHARED SchroedingerGreen1D_lib.cpp )
add_library(FitREQ_lib               SHARED FitREQ_lib.cpp )

