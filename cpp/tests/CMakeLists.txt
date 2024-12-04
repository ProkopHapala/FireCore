include(FetchContent)

FetchContent_Declare(
  googletest
  GIT_REPOSITORY https://github.com/google/googletest.git
  GIT_TAG        release-1.11.0
)
FetchContent_MakeAvailable(googletest)

# Add the executable for the test
add_executable(test_add test_add.cpp)

# Link against gtest and gtest_main directly
target_link_libraries(test_add PRIVATE 
    gtest 
    gtest_main
)

# Register the test with CTest
add_test(NAME AddTest COMMAND "${CMAKE_BINARY_DIR}/tests/test_add")


