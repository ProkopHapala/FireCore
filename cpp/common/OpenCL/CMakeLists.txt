
#add_library( MoleculeWorld2D OBJECT MoleculeWorld2D.cpp )
#add_library( MolecularEngine OBJECT MolecularWorld.cpp MoleculeType.cpp AtomTypes.cpp  )

file(GLOB OCL_HEADERS_ *.h )
add_custom_target(OCL_HEADERS SOURCES ${OCL_HEADERS_} )
