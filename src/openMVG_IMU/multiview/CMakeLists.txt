
file(
  GLOB_RECURSE
  multiview_files_header
  *.hpp
)
file(
  GLOB_RECURSE
  multiview_files_cpp
  *.cpp
)


#Remove the unit test files (not been used by the library)


add_library(openMVG_IMU_multiview
  ${multiview_files_header}
  ${multiview_files_cpp})
target_compile_features(openMVG_IMU_multiview INTERFACE ${CXX11_FEATURES})
set_property(TARGET openMVG_IMU_multiview PROPERTY FOLDER OpenMVG/OpenMVG_IMU)
install(TARGETS openMVG_IMU_multiview DESTINATION lib EXPORT openMVG-targets)
target_link_libraries(openMVG_IMU_multiview
  PUBLIC
    openMVG_multiview
    openMVG_numeric
  PRIVATE
    openMVG_graph
    ${CERES_LIBRARIES}
)
target_include_directories(openMVG_IMU_multiview
  PRIVATE
    ${CERES_INCLUDE_DIRS}
  PUBLIC
    $<INSTALL_INTERFACE:include>
)
set_target_properties(openMVG_IMU_multiview PROPERTIES SOVERSION ${OPENMVG_VERSION_MAJOR} VERSION "${OPENMVG_VERSION_MAJOR}.${OPENMVG_VERSION_MINOR}")



