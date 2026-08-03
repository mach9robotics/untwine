
set(GOOGLETEST_VERSION 1.12.1)

# Since tests link CRT dynamically (/MD[d]), require gtest to
#link dynamically too (default is /MT[d])
option(gtest_force_shared_crt "Always use shared Visual C++ run-time DLL" ON)

add_subdirectory(gtest)

###############################################################################
# Add a test target.
# _name The driver name.
# ARGN :
#    FILES the source files for the test
#    LINK_WITH link test executable with libraries
#    INCLUDES header file directories
#
macro(UNTWINE_ADD_TEST _name)
    set(options)
    set(oneValueArgs)
    set(multiValueArgs FILES LINK_WITH INCLUDES SYSTEM_INCLUDES)
    cmake_parse_arguments(UNTWINE_ADD_TEST
        "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})
    add_executable(${_name}
        ${UNTWINE_ADD_TEST_FILES}
        ${LAZPERF_SRCS}
    )
    untwine_target_compile_settings(${_name})
    # gtest/include MUST precede PDAL's include dir: some PDAL installs (e.g.
    # conda) ship their own <gtest/gtest.h> of a different version, and picking
    # that header up while linking the bundled libgtest.a yields undefined-symbol
    # link errors from the mismatched MakeAndRegisterTestInfo signature.
    target_include_directories(${_name} PRIVATE
        ${PROJECT_SOURCE_DIR}
        ${PROJECT_SOURCE_DIR}/${UNTWINE_OS_DIR}
        ${CMAKE_CURRENT_SOURCE_DIR}/gtest/include
        ${PDAL_INCLUDE_DIRS}
        ${CMAKE_CURRENT_BINARY_DIR}
    )
    set_property(TARGET ${_name} PROPERTY FOLDER "Tests")
    target_link_libraries(${_name}
        PRIVATE
            gtest
            ${PDAL_LIBRARIES}
            ${UNTWINE_ADD_TEST_LINK_WITH}
    )
    add_test(NAME ${_name}
        COMMAND
            "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${_name}"
        WORKING_DIRECTORY
            "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/.."
    )
endmacro()

configure_file(TestConfig.hpp.in ${CMAKE_CURRENT_BINARY_DIR}/TestConfig.hpp)

UNTWINE_ADD_TEST(untwine_test FILES Tests.cpp)

# Pure-logic unit tests for each pipeline step + the seams between them. These
# compile the algorithmic source units directly (Grid for prep/EPF binning,
# Stats for the BU rollup); the threaded orchestration classes are exercised
# end-to-end by untwine_test above.
UNTWINE_ADD_TEST(untwine_unit_test
    FILES
        UnitTestMain.cpp
        VoxelKeyTest.cpp
        GridKeyTest.cpp
        GridTest.cpp
        PointTest.cpp
        StatsTest.cpp
        CommonTest.cpp
        PipelineIntegrationTest.cpp
        ChunkBuildTest.cpp
        ChunkLutTest.cpp
        ReChunkTest.cpp
        ThreadPoolTest.cpp
        ${PROJECT_SOURCE_DIR}/pointio/Grid.cpp
        ${PROJECT_SOURCE_DIR}/bu/Stats.cpp
        ${PROJECT_SOURCE_DIR}/bu/ChunkPlan.cpp
        ${PROJECT_SOURCE_DIR}/chunker/ChunkLut.cpp
        ${PROJECT_SOURCE_DIR}/chunker/ReChunk.cpp
        ${PROJECT_SOURCE_DIR}/untwine/ThreadPool.cpp
)
