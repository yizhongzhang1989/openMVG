
file(
  GLOB_RECURSE
  utils_files_header
  *.hpp
)
file(
  GLOB_RECURSE
  utils_files_cpp
  *.cpp
)


add_library(openMVG_IMU_utils ${utils_files_header} ${utils_files_cpp})
target_compile_features(openMVG_IMU_utils INTERFACE ${CXX11_FEATURES})
target_link_libraries(openMVG_IMU_utils
  PUBLIC
    openMVG_sfm
    
    ${OPENMVG_LIBRARY_DEPENDENCIES}
)

target_link_libraries(openMVG_IMU_utils
  LINK_PRIVATE
    
    ${STLPLUS_LIBRARY}
)
target_include_directories(openMVG_IMU_utils
  PRIVATE
    
  PUBLIC
    $<INSTALL_INTERFACE:include>
    $<INSTALL_INTERFACE:include/openMVG>
)
if (MSVC)
  set_target_properties(openMVG_IMU_utils PROPERTIES COMPILE_FLAGS "/bigobj")
endif (MSVC)
set_target_properties(openMVG_IMU_utils PROPERTIES SOVERSION ${OPENMVG_VERSION_MAJOR} VERSION "${OPENMVG_VERSION_MAJOR}.${OPENMVG_VERSION_MINOR}")
set_property(TARGET openMVG_IMU_utils PROPERTY FOLDER OpenMVG/OpenMVG_IMU)
install(TARGETS openMVG_IMU_utils DESTINATION lib EXPORT openMVG-targets)
