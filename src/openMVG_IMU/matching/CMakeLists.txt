
project(matching)

file(
  GLOB
  matching_files_header
  *.hpp
)
file(
  GLOB
  matching_files_cpp
  *.cpp
)


set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)

set_source_files_properties(${matching_files_cpp} PROPERTIES LANGUAGE CXX)
add_library(openMVG_IMU_matching
  ${matching_files_header}
  ${matching_files_cpp})
target_compile_features(openMVG_IMU_matching INTERFACE ${CXX11_FEATURES})
target_include_directories(openMVG_IMU_matching
  PUBLIC
    $<BUILD_INTERFACE:${FLANN_INCLUDE_DIRS}>
    $<INSTALL_INTERFACE:include>
    $<INSTALL_INTERFACE:include/openMVG>
)
target_link_libraries(openMVG_IMU_matching
  PRIVATE
    ${STLPLUS_LIBRARY}
    ${OPENMVG_LIBRARY_DEPENDENCIES}
  PUBLIC
    openMVG_features
    openMVG_matching
    Threads::Threads
    cereal
)
if (NOT DEFINED OpenMVG_USE_INTERNAL_FLANN)
target_link_libraries(openMVG_IMU_matching
  PUBLIC
    ${FLANN_LIBRARIES}
)
endif()
set_target_properties(openMVG_IMU_matching PROPERTIES SOVERSION ${OPENMVG_VERSION_MAJOR} VERSION "${OPENMVG_VERSION_MAJOR}.${OPENMVG_VERSION_MINOR}")
set_property(TARGET openMVG_IMU_matching PROPERTY FOLDER OpenMVG/OpenMVG_IMU)

if (USE_AVX2)
  target_compile_options(openMVG_IMU_matching PUBLIC "-DOPENMVG_USE_AVX2")
  if (UNIX)
    target_compile_options(openMVG_IMU_matching PUBLIC "-mavx2")
  endif (UNIX)
endif (USE_AVX2)

install(TARGETS openMVG_IMU_matching DESTINATION lib EXPORT openMVG-targets)
