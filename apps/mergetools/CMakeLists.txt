# apps/mergetools/CMakeLists.txt


add_executable(mergeRinNav mergeRinNav.cpp)


add_executable(mergeRinMet mergeRinMet.cpp)

if (CMAKE_SYSTEM_NAME MATCHES "SunOS")
  target_link_libraries(mergeRinMet gen)
  target_link_libraries(mergeRinNav gen)
endif()

