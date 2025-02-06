# Set the output directory explicitly
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/results)

# Add the executable for the test
add_executable(free_energy free_energy.cpp)
target_link_libraries(free_energy PRIVATE MMFF_lib)
