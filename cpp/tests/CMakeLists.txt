include(FetchContent)

FetchContent_Declare(
  googletest
  GIT_REPOSITORY https://github.com/google/googletest.git
  GIT_TAG        release-1.11.0
)
FetchContent_MakeAvailable(googletest)
# Set the output directory explicitly
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/tests)

# Add the executable for the test
add_executable(test_add test_add.cpp)
target_link_libraries(test_add PRIVATE 
   gtest 
   gtest_main
)

add_executable(test_free_energy test_free_energy.cpp)
target_link_libraries(test_free_energy PRIVATE MMFF_lib gtest gtest_main)


add_executable(test_MMFFsp3_loc test_MMFFsp3_loc.cpp)
target_link_libraries(test_MMFFsp3_loc PRIVATE MMFF_lib gtest gtest_main)


# Register the test with CTest using the correct path
add_test(NAME AddTest COMMAND $<TARGET_FILE:test_add>)
add_test(NAME FreeEnergyTest COMMAND $<TARGET_FILE:test_free_energy>)
add_test(NAME MMFFsp3_locTest COMMAND $<TARGET_FILE:test_MMFFsp3_loc>)
