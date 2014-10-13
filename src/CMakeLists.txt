CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

SET(FTMAKE_MINIMUM_REQUIRED_VERSION 1.7)

PROJECT(CH4)
set (CH4_VERSION_MAJOR 1)
set (CH4_VERSION_MINOR 0)

ENABLE_TESTING()

set(TARGET_NAME ch4)

set(DOC_DIRECTORY doc)

execute_process(
        COMMAND git rev-parse --verify --symbolic-full-name --abbrev-ref HEAD
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        OUTPUT_VARIABLE git-branch
        OUTPUT_STRIP_TRAILING_WHITESPACE
)

message(STATUS "FTMAKE: ${git-branch}")

if (ON_MASTER_TEST_FTMAKE_DEVELOP)
message(STATUS "FTMAKE: Nightly-ftmake-rc mode is on")
set(FTMAKE "ftmake-develop")
endif()

# Try to find the ftmake directories
if(FTMAKE)
else(FTMAKE)
set(FTMAKE "ftmake")
endif(FTMAKE)

message(STATUS "FTMAKE: ${FTMAKE}")


SET(TRIAL_PATHS
  $ENV{HOME}/${FTMAKE}
  $ENV{HOME}/projects/${FTMAKE}
  /mnt/ernie3/reijers/${FTMAKE}
  NO_DEFAULT_PATH
)



FIND_PATH(FTMAKE_SOURCE_DIR CMakeLists.txt ${TRIAL_PATHS})


# Launch cmake
if(FTMAKE_SOURCE_DIR)
include(${FTMAKE_SOURCE_DIR}/CMakeLists.txt)
ELSE(FTMAKE_SOURCE_DIR)
message(STATUS "${FTMAKE_SOURCE_DIR}/CMakeLists.txt not found")
message(STATUS "Go in your home and then:")
message(STATUS "git clone gitolite@bert.phys.tue.nl:toschi/ftmake")
message(STATUS "after that you can cmake . in this directory")
message(FATAL_ERROR "aborting")
ENDIF(FTMAKE_SOURCE_DIR)


message(STATUS "${git-branch}")

set(BUILDNAME "${git-branch}@${git-short-hash} (${git-ftmake-branch}@${git-ftmake-short-hash})")

INCLUDE(CTest)
message("CTEST_FULL_OUTPUT")




