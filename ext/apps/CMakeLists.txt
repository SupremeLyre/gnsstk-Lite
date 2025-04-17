# apps/CMakeLists.txt
link_libraries(gnsstk)
include_directories(${PROJECT_SOURCE_DIR}/include)
add_subdirectory (filetools)
add_subdirectory (geomatics)

# If on UNIX, do these
if (UNIX)
    add_subdirectory (rfw)
endif (UNIX)


