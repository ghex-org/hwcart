find_package(MPI REQUIRED)

set_target_properties(hwcart PROPERTIES PUBLIC_HEADER "${PROJECT_SOURCE_DIR}/include/hwcart/hwcart.h;${PROJECT_SOURCE_DIR}/include/hwcart/hwcart_utils.h")

if (USE_HWLOC)
   message(STATUS "linking HWLOC " ${HWLOC_INCLUDE_DIRS})
   find_package(HWLOC REQUIRED)
   target_compile_definitions(hwcart PUBLIC USE_HWLOC)
   target_link_libraries (hwcart PUBLIC HWLOC::libhwloc)
   target_sources(hwcart PUBLIC hwcart.c)
else ()
   target_sources(hwcart PUBLIC hwcart_ompi.c)
endif (USE_HWLOC)

target_sources(hwcart PUBLIC hwcart_common.c)

target_include_directories(
  hwcart
  PUBLIC
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include/>
  $<INSTALL_INTERFACE:include>
)

install (TARGETS hwcart 
	 LIBRARY DESTINATION lib
	 PUBLIC_HEADER DESTINATION include)
