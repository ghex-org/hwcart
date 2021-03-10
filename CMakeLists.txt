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

find_package(MPI REQUIRED)

# do we use HWLOC?
option (HWCART_USE_HWLOC "use hwloc" OFF)
if (HWCART_USE_HWLOC)
   find_package(HWLOC)
   if (HWLOC_FOUND)
     list (APPEND DEPENDENCY_TARGETS HWLOC::libhwloc)
     message(STATUS "Using HWLOC hardware locality features")
   endif (HWLOC_FOUND)
else ()
   message(STATUS "Using OpenMPI hardware locality features. Can only compile using OpenMPI.")
endif (HWCART_USE_HWLOC)

# main shared lib
add_library(hwcart SHARED)
add_library(HWCART::hwcart ALIAS hwcart)
target_include_directories(hwcart PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>)
target_link_libraries (hwcart PUBLIC MPI::MPI_C)

# Fortran bindings
option (HWCART_USE_FORTRAN "Enable Fortran bindings" OFF)
if (HWCART_USE_FORTRAN)
  enable_language(Fortran)
  add_library(hwcart_f SHARED)
  target_link_libraries (hwcart_f PUBLIC hwcart)
endif(HWCART_USE_FORTRAN)

add_subdirectory(src)

option(HWCART_USE_TESTS "build tests" OFF)
if (HWCART_USE_TESTS)
  add_subdirectory(tests)
endif()

include(GNUInstallDirs)
install(TARGETS hwcart
    EXPORT HWCART-targets
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
install(DIRECTORY include/
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
install(EXPORT HWCART-targets
  FILE HWCART-targets.cmake
  NAMESPACE HWCART::
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake)
include(CMakePackageConfigHelpers)
configure_package_config_file(${CMAKE_CURRENT_SOURCE_DIR}/cmake/HWCARTConfig.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/HWCARTConfig.cmake
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake)
write_basic_package_version_file(HWCARTConfigVersion.cmake
    VERSION ${PROJECT_VERSION} COMPATIBILITY SameMajorVersion)
install(
    FILES
        ${CMAKE_CURRENT_BINARY_DIR}/HWCARTConfig.cmake
        ${CMAKE_CURRENT_BINARY_DIR}/HWCARTConfigVersion.cmake
        ${CMAKE_CURRENT_LIST_DIR}/cmake/FindHWLOC.cmake
    DESTINATION
        ${CMAKE_INSTALL_LIBDIR}/cmake)
export(EXPORT HWCART-targets
  FILE "${CMAKE_CURRENT_BINARY_DIR}/HWCART-targets.cmake")
