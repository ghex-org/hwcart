set(_tests test_hwcart test_hwcart2mpicart)
foreach (_t ${_tests})
  add_executable(${_t} ${_t}.c)

  # mark RPATH to be used when installing, but only if it's not a system directory
  list(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
  if("${isSystemDir}" STREQUAL "-1")
    set_property(TARGET ${_t} PROPERTY INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
  endif("${isSystemDir}" STREQUAL "-1")

  target_link_libraries (${_t} PUBLIC ${DEPENDENCY_TARGETS})
  target_link_libraries (${_t} PUBLIC MPI::MPI_C)
  target_link_libraries (${_t} PUBLIC hwcart)

  install (TARGETS ${_t} DESTINATION bin)
endforeach()

if (HWCART_USE_FORTRAN)
  foreach (_t ${_tests})
    add_executable(${_t}_f ${_t}.f90 )
    target_include_directories(${_t}_f PUBLIC
      $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}>
      $<INSTALL_INTERFACE:include>)

    # mark RPATH to be used when installing, but only if it's not a system directory
    list(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
    if("${isSystemDir}" STREQUAL "-1")
      set_property(TARGET ${_t} PROPERTY INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
    endif("${isSystemDir}" STREQUAL "-1")

    target_link_libraries (${_t}_f PUBLIC ${DEPENDENCY_TARGETS})
    target_link_libraries (${_t}_f PUBLIC MPI::MPI_Fortran)
    target_link_libraries (${_t}_f PUBLIC hwcart)

    install (TARGETS ${_t}_f DESTINATION bin)
  endforeach()
endif (HWCART_USE_FORTRAN)  

