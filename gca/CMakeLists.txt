cmake_minimum_required(VERSION 3.10)
project(GCATreeProject CXX)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Add debugging symbols
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -g")

# --- Include Directories ---
# All headers are in the include directory
include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)

# --- Find Dependencies (PBC, GMP, and GTest) ---

find_package(PBC REQUIRED)
find_package(GMP REQUIRED)
# find_package(GTest REQUIRED)  # Add GTest

# Include directories for dependencies
include_directories(${PBC_INCLUDE_DIRS})
include_directories(${GMP_INCLUDE_DIRS})

# --- Source Files ---
# Get all .cpp files from subdirectories in src/
file(GLOB SUBDIR_SRCS "src/*/*.cpp")
# Also include any .cpp files directly in the src directory
file(GLOB SRC_ROOT "src/*.cpp")

set(ALL_SRCS
    ${SUBDIR_SRCS}
    ${SRC_ROOT}
)

# --- Merkle Tree Library ---
# Create a library for Merkle tree functionality
add_library(mht_lib STATIC
    src/mht/merkle.cpp
    src/mht/hash.cpp
    src/mht/transaction.cpp
)

# --- Accumulator Library ---
add_library(acc_lib STATIC
    src/acc/acc_value.cpp
    src/acc/accumulator.cpp
    src/acc/set.cpp
    src/acc/keys.cpp
    src/acc/keys_io.cpp
    src/acc/proof.cpp
)

# --- GCA Library ---
add_library(gca_lib STATIC
    src/gca/gca_tree.cpp
    src/gca/sorted_list.cpp
    src/gca/prefix_set.cpp
    src/gca/triplet.cpp
)

# Set include directories for the library
target_include_directories(mht_lib PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)

target_include_directories(acc_lib PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${PBC_INCLUDE_DIRS}
    ${GMP_INCLUDE_DIRS}
)

target_include_directories(gca_lib PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${PBC_INCLUDE_DIRS}
    ${GMP_INCLUDE_DIRS}
)

# Link PBC and GMP libraries to acc_lib
target_link_libraries(acc_lib PUBLIC
    ${PBC_LIBRARIES}
    ${GMP_LIBRARIES}
)

# Link dependencies to gca_lib
target_link_libraries(gca_lib PUBLIC
    mht_lib
    acc_lib
    ${PBC_LIBRARIES}
    ${GMP_LIBRARIES}
)

# --- Executable Target ---
add_executable(gca_tree ${ALL_SRCS})

# --- Link Libraries ---
target_link_libraries(gca_tree
    PRIVATE
    mht_lib  # Add Merkle tree library
    acc_lib  # Add Accumulator library
    gca_lib  # Add GCA library
    ${PBC_LIBRARIES}  # Move PBC libraries to the end
    ${GMP_LIBRARIES}
)

# --- Enable Testing ---
enable_testing()

# # --- Configure CTest for debugging ---
# set(CTEST_DEBUG_OPTIONS "--debug-output --output-log=ctest_debug.log")
# set(CTEST_DEBUG_OPTIONS "${CTEST_DEBUG_OPTIONS} --debug-output --output-on-failure")
# set(CTEST_DEBUG_OPTIONS "${CTEST_DEBUG_OPTIONS} --debug-output --verbose")
# set(CTEST_DEBUG_OPTIONS "${CTEST_DEBUG_OPTIONS} --debug-output --stop-on-failure")

# # Set CTest to use lldb for debugging
# set(CTEST_DEBUG_OPTIONS "${CTEST_DEBUG_OPTIONS} --debug-output --debugger=lldb")

# # --- Add Test Subdirectory ---
# add_subdirectory(test)

# --- Output ---
message(STATUS "Found PBC: ${PBC_FOUND}")
message(STATUS "PBC Include Dirs: ${PBC_INCLUDE_DIRS}")
message(STATUS "PBC Libraries: ${PBC_LIBRARIES}")
message(STATUS "Found GMP: ${GMP_FOUND}")
message(STATUS "GMP Include Dirs: ${GMP_INCLUDE_DIRS}")
message(STATUS "GMP Libraries: ${GMP_LIBRARIES}")
message(STATUS "Found GTest: ${GTEST_FOUND}")

# set(CTEST_DISABLE_TESTS TRUE)

# # Add a build target that doesn't run tests
# add_custom_target(build COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR})
