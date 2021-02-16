cmake_minimum_required(VERSION 3.14.5)

file(STRINGS "version.txt" _hwcart_version)
project(HWCart VERSION ${_hwcart_version} LANGUAGES C)
unset(_hwcart_version)

if(NOT CMAKE_BUILD_TYPE) # AND NOT CMAKE_CONFIGURATION_TYPES)
   set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build." FORCE)
   # Set the possible values of build type for cmake-gui
   set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release")
   set_property(CACHE CMAKE_C_FLAGS PROPERTY VALUE "-Wall -Wextra -Wpedantic -Wno-unknown-pragmas")
   set_property(CACHE CMAKE_C_FLAGS_RELEASE PROPERTY VALUE "-O3 -DNDEBUG -Wno-unused-parameter")
endif()

add_compile_definitions(_GNU_SOURCE)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")

# do we use HWLOC?
option (USE_HWLOC "use hwloc" OFF)
if (USE_HWLOC)
   find_package(HWLOC)
   if (HWLOC_FOUND)
     list (APPEND DEPENDENCY_TARGETS HWLOC::libhwloc)
     message(STATUS "Using HWLOC hardware locality features")
   endif (HWLOC_FOUND)
else ()
   message(STATUS "Using OpenMPI hardware locality features. Can only compile using OpenMPI.")
endif (USE_HWLOC)

add_library(hwcart SHARED)

add_subdirectory(src)

option (USE_TESTS "build tests" OFF)
if (USE_TESTS)
  add_subdirectory(tests)
endif()
