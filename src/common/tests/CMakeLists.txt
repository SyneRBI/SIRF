add_executable(test_version test_version.cpp)
target_include_directories(test_version PUBLIC ${CMAKE_BINARY_DIR}/cmake/include)

add_test(NAME SIRF_TEST_COMMON_CPLUSPLUS
         COMMAND test_version )