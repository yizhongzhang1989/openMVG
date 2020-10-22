#
#
#set(OpenMVG_USE_INTERNAL_CERES ON)
#set(CERES_INCLUDE_DIRS
#        ${CMAKE_CURRENT_SOURCE_DIR}/third_party/ceres-solver/include
#        ${CMAKE_CURRENT_SOURCE_DIR}/third_party/ceres-solver/internal/ceres/miniglog
#        ${PROJECT_BINARY_DIR}/third_party/ceres-solver/config)
#
#STRING(REGEX REPLACE "version ([0-9.]+).*" "\\1" CERES_VERSION ${CERES_CONFIG})
#set(CERES_LIBRARIES openMVG_ceres)
#
#include_directories(${CERES_INCLUDE_DIRS})

file(
        GLOB_RECURSE
        vi_sfm_files_header
        *.hpp
)
file(
        GLOB_RECURSE
        vi_sfm_files_cpp
        *.cpp
)
#file(GLOB_RECURSE REMOVEFILESUNITTEST *_test.cpp)
#Remove the future main files
#list(REMOVE_ITEM sfm_files_cpp ${REMOVEFILESUNITTEST})
#Remove the test header main files
#file(GLOB_RECURSE REMOVEFILESUNITTEST *_test.hpp)
#list(REMOVE_ITEM sfm_files_header ${REMOVEFILESUNITTEST})


add_library(openMVG_vi_sfm ${vi_sfm_files_header} ${vi_sfm_files_cpp})
target_compile_features(openMVG_vi_sfm INTERFACE ${CXX11_FEATURES})
target_link_libraries(openMVG_vi_sfm
        PUBLIC
        openMVG_geometry
        openMVG_features
        openMVG_graph
        openMVG_matching
        openMVG_multiview
        cereal
        ${OPENMVG_LIBRARY_DEPENDENCIES}
        )

target_link_libraries(openMVG_vi_sfm
        LINK_PRIVATE
        ${CERES_LIBRARIES}
        openMVG_graph
        openMVG_image
        openMVG_lInftyComputerVision
        openMVG_system
        openMVG_sfm
        ${STLPLUS_LIBRARY}
        )
target_include_directories(openMVG_vi_sfm
        PRIVATE
        ${CERES_INCLUDE_DIRS}
        PUBLIC
        $<INSTALL_INTERFACE:include>
        $<INSTALL_INTERFACE:include/openMVG>
        )
if (MSVC)
    set_target_properties(openMVG_vi_sfm PROPERTIES COMPILE_FLAGS "/bigobj")
endif (MSVC)
set_target_properties(openMVG_vi_sfm PROPERTIES SOVERSION ${OPENMVG_VERSION_MAJOR} VERSION "${OPENMVG_VERSION_MAJOR}.${OPENMVG_VERSION_MINOR}")
set_property(TARGET openMVG_vi_sfm PROPERTY FOLDER OpenMVG/OpenMVG)
install(TARGETS openMVG_vi_sfm DESTINATION lib EXPORT openMVG-targets)

UNIT_TEST(openMVG sfm_data_io "openMVG_vi_sfm;${STLPLUS_LIBRARY}")
UNIT_TEST(openMVG sfm_data_BA "openMVG_multiview_test_data;openMVG_vi_sfm;${STLPLUS_LIBRARY}")
UNIT_TEST(openMVG sfm_data_utils "openMVG_vi_sfm;${STLPLUS_LIBRARY}")
UNIT_TEST(openMVG sfm_data_filters "openMVG_vi_sfm")
UNIT_TEST(openMVG sfm_data_graph_utils "openMVG_vi_sfm")
UNIT_TEST(openMVG sfm_data_triangulation "openMVG_vi_sfm;openMVG_multiview_test_data")

