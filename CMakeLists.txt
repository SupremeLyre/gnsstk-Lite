# ============================================================
# Name    = $GNSSTK/CMakeLists.txt
# Purpose = Generator of build framework (e.g. Makefiles) for GNSSTk
# Notes   = This is the top-level CMake input file
# Depends on $GNSSTK/BuildSetup.cmake
# Is dependend on by $GNSSTK/build.sh
# ============================================================

cmake_minimum_required(VERSION 3.7.2)

project(GNSSTK)
set(GNSSTK_VERSION_MAJOR "14")
set(GNSSTK_VERSION_MINOR "3")
set(GNSSTK_VERSION_PATCH "0")
set(GNSSTK_VERSION "${GNSSTK_VERSION_MAJOR}.${GNSSTK_VERSION_MINOR}.${GNSSTK_VERSION_PATCH}")

# Test data file source directory
set(GNSSTK_TEST_DATA_DIR ${PROJECT_SOURCE_DIR}/data)

# Test output directory
set(GNSSTK_TEST_OUTPUT_DIR ${PROJECT_BINARY_DIR}/Testing/Temporary)

# This sets up variables contining GNU standard installation locations.
include(GNUInstallDirs)

# Set a filename for collecting exported targets.
set(EXPORT_TARGETS_FILENAME "GNSSTKTargets")

if(WIN32)
    # Make cmake do most of the work of creating the DLL exports
    set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS TRUE)

    # put everything in the same directory to make it easier for tests
    # to find the DLL.
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
else(WIN32)
    # This disables the Windows-specific export things for non-Windows platforms.
    add_definitions(-DGNSSTK_STATIC_DEFINE)
endif(WIN32)

include(GenerateExportHeader)

# ============================================================
# Core Library Target Files
# ============================================================

# ----------------------------------------
# Define $GNSSTK/core/ library source files
# ----------------------------------------
file(GLOB_RECURSE CORE_SRC_FILES "core/lib/*.cpp" "core/lib/*.c")
file(GLOB_RECURSE CORE_INC_FILES "core/lib/*.h" "core/lib/*.hpp")

# ----------------------------------------
# Define /core library include directories
# ----------------------------------------

# initialize list of include directories
set(CORE_INC_DIRS "")

foreach(_headerFile ${CORE_INC_FILES})
    get_filename_component(_dir ${_headerFile} PATH)
    list(APPEND CORE_INC_DIRS ${_dir})
endforeach()

list(REMOVE_DUPLICATES CORE_INC_DIRS)

# Add every directory containing a header file
# to the project(gnsstk) include_directories
include_directories(${CORE_INC_DIRS})

# define src and include files needed to build library target
set(GNSSTK_SRC_FILES "")
set(GNSSTK_INC_FILES "")
list(APPEND GNSSTK_SRC_FILES ${CORE_SRC_FILES})
list(APPEND GNSSTK_INC_FILES ${CORE_INC_FILES})

# Remove getopt.h from non-Windows installs
if(NOT WIN32)
    foreach(_headerFile ${GNSSTK_INC_FILES})
        get_filename_component(_name ${_headerFile} NAME)

        if(${_name} MATCHES "getopt.h")
            list(REMOVE_ITEM GNSSTK_INC_FILES ${_headerFile})
        endif()
    endforeach()
endif()

# ============================================================
# Define $GNSSTK/ext/ additions to Library Target Files
# ============================================================
add_definitions(-DBUILD_EXT)
file(GLOB_RECURSE EXT_SRC_FILES "ext/lib/*.cpp" "ext/lib/*.c")
file(GLOB_RECURSE EXT_INC_FILES "ext/lib/*.h" "ext/lib/*.hpp")

# Define ext library include directories
set(EXT_INC_DIRS "")

foreach(_headerFile ${EXT_INC_FILES})
    get_filename_component(_dir ${_headerFile} PATH)
    get_filename_component(_name ${_headerFile} NAME)
    list(APPEND EXT_INC_DIRS ${_dir})
endforeach()

list(REMOVE_DUPLICATES EXT_INC_DIRS)

# Add every directory containing a header file
# to the project(gnsstk) include_directories
include_directories(${EXT_INC_DIRS})

# append ext include files needed to build library target
list(APPEND GNSSTK_SRC_FILES ${EXT_SRC_FILES})
list(APPEND GNSSTK_INC_FILES ${EXT_INC_FILES})

# ============================================================
# GNSSTK Library, Build and Install Targets
# ============================================================
if(${PROFILER})
    # don't build the shared library with profiler enabled
    add_library(gnsstk STATIC ${GNSSTK_SRC_FILES} ${GNSSTK_INC_FILES})
elseif(WIN32)
    add_library(gnsstk ${GNSSTK_SRC_FILES} ${GNSSTK_INC_FILES})
else()
    add_library(gnsstk SHARED ${GNSSTK_SRC_FILES} ${GNSSTK_INC_FILES})
endif()

# always generate the header because it's an include file whose
# absence would break the build on non-windows.
generate_export_header(gnsstk)
include_directories(AFTER SYSTEM ${PROJECT_BINARY_DIR})
