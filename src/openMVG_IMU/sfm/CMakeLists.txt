
file(
  GLOB_RECURSE
  sfm_files_header
  *.hpp
)
file(
  GLOB_RECURSE
  sfm_files_cpp
  *.cpp
)


add_library(openMVG_IMU_sfm ${sfm_files_header} ${sfm_files_cpp})
target_compile_features(openMVG_IMU_sfm INTERFACE ${CXX11_FEATURES})
target_link_libraries(openMVG_IMU_sfm
  PUBLIC
    openMVG_sfm
    openMVG_IMU_multiview
    openMVG_IMU_matching_image_collection
    openMVG_geometry
    openMVG_features
    openMVG_graph
    openMVG_matching
    openMVG_multiview
    cereal
    ${OPENMVG_LIBRARY_DEPENDENCIES}
)

target_link_libraries(openMVG_IMU_sfm
  LINK_PRIVATE
    ${CERES_LIBRARIES}
    openMVG_graph
    openMVG_image
    openMVG_lInftyComputerVision
    openMVG_system
    ${STLPLUS_LIBRARY}
)
target_include_directories(openMVG_IMU_sfm
  PRIVATE
    ${CERES_INCLUDE_DIRS}
  PUBLIC
    $<INSTALL_INTERFACE:include>
    $<INSTALL_INTERFACE:include/openMVG>
)
if (MSVC)
  set_target_properties(openMVG_IMU_sfm PROPERTIES COMPILE_FLAGS "/bigobj")
endif (MSVC)
set_target_properties(openMVG_IMU_sfm PROPERTIES SOVERSION ${OPENMVG_VERSION_MAJOR} VERSION "${OPENMVG_VERSION_MAJOR}.${OPENMVG_VERSION_MINOR}")
set_property(TARGET openMVG_IMU_sfm PROPERTY FOLDER OpenMVG/OpenMVG_IMU)
install(TARGETS openMVG_IMU_sfm DESTINATION lib EXPORT openMVG-targets)

UNIT_TEST(openMVG sfm_data_io "openMVG_IMU_sfm;${STLPLUS_LIBRARY}")
UNIT_TEST(openMVG sfm_data_BA "openMVG_multiview_test_data;openMVG_IMU_sfm;${STLPLUS_LIBRARY}")
UNIT_TEST(openMVG sfm_data_utils "openMVG_IMU_sfm;${STLPLUS_LIBRARY}")
UNIT_TEST(openMVG sfm_data_filters "openMVG_IMU_sfm")
UNIT_TEST(openMVG sfm_data_graph_utils "openMVG_IMU_sfm")
UNIT_TEST(openMVG sfm_data_triangulation "openMVG_IMU_sfm;openMVG_multiview_test_data")

add_subdirectory(pipelines)
