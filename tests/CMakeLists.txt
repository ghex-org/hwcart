add_executable(hwcart_test test.c)

# mark RPATH to be used when installing, but only if it's not a system directory
list(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
if("${isSystemDir}" STREQUAL "-1")
  set_property(TARGET hwcart_test PROPERTY INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
endif("${isSystemDir}" STREQUAL "-1")

target_link_libraries (hwcart_test PUBLIC ${DEPENDENCY_TARGETS})
target_link_libraries (hwcart_test PUBLIC MPI::MPI_C)
target_link_libraries (hwcart_test PUBLIC hwcart)

install (TARGETS hwcart_test DESTINATION bin)
