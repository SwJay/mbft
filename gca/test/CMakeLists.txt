cmake_minimum_required(VERSION 3.10)
project(AccumulatorTest)

# Find required packages
find_package(GTest REQUIRED)
find_package(PBC REQUIRED)
find_package(GMP REQUIRED)

# Include directories - adjust if your mht headers are elsewhere
# Assumes mht headers are accessible via 'mht/...' relative to the project root
include_directories(${CMAKE_SOURCE_DIR}/include) # Or specify the exact path to the dir containing 'mht'
include_directories(${PBC_INCLUDE_DIRS})
include_directories(${GMP_INCLUDE_DIRS})

# Add the test executable
add_executable(merkle_test merkle_test.cpp)
add_executable(nested_difference_test nested_difference_test.cpp)
add_executable(test_union test_union.cpp)
add_executable(test_nested_union test_nested_union.cpp)
add_executable(test_compare_keys test_compare_keys.cpp)

# Link the test executable to GTest and your Merkle Tree library
target_link_libraries(merkle_test PRIVATE GTest::GTest GTest::Main mht_lib)
target_link_libraries(nested_difference_test PRIVATE 
    GTest::GTest 
    GTest::Main 
    acc_lib 
    mht_lib
    ${PBC_LIBRARIES}
    ${GMP_LIBRARIES}
)
target_link_libraries(test_union PRIVATE
    GTest::GTest
    GTest::Main
    acc_lib
    mht_lib
    ${PBC_LIBRARIES}
    ${GMP_LIBRARIES}
)
target_link_libraries(test_nested_union PRIVATE
    GTest::GTest
    GTest::Main
    acc_lib
    mht_lib
    ${PBC_LIBRARIES}
    ${GMP_LIBRARIES}
)
target_link_libraries(test_compare_keys PRIVATE
    acc_lib
    mht_lib
    ${PBC_LIBRARIES}
    ${GMP_LIBRARIES}
)

# Enable testing
include(GoogleTest)
gtest_discover_tests(merkle_test)
gtest_discover_tests(nested_difference_test)
gtest_discover_tests(test_union)
gtest_discover_tests(test_nested_union)
 