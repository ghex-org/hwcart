add_executable(test test.c)

find_package(MPI REQUIRED)

# mark RPATH to be used when installing, but only if it's not a system directory
list(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
if("${isSystemDir}" STREQUAL "-1")
  set_property(TARGET test PROPERTY INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
endif("${isSystemDir}" STREQUAL "-1")

target_link_libraries (test PUBLIC ${DEPENDENCY_TARGETS})
target_link_libraries (test PUBLIC MPI::MPI_C)
target_link_libraries (test PUBLIC hwcart)

install (TARGETS test DESTINATION bin)
