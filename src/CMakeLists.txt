if (HWCART_USE_HWLOC)
   message(STATUS "linking HWLOC " ${HWLOC_INCLUDE_DIRS})
   find_package(HWLOC REQUIRED)
   target_compile_definitions(hwcart PUBLIC HWCART_USE_HWLOC)
   target_link_libraries (hwcart PUBLIC HWLOC::libhwloc)
   target_sources(hwcart PRIVATE hwcart.c)
else ()
   target_sources(hwcart PRIVATE hwcart_ompi.c)
endif (HWCART_USE_HWLOC)

if (HWCART_USE_FORTRAN)
  target_sources(hwcart PRIVATE hwcart_mod.f90)
  target_include_directories (hwcart PUBLIC ${MPI_Fortran_INCLUDE_DIRS})
  target_compile_definitions(hwcart PUBLIC HWCART_USE_FORTRAN)
endif(HWCART_USE_FORTRAN)

target_sources(hwcart PRIVATE hwcart_common.c)
