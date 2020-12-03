add_executable(openMVG_main_AppendIncrementalSfM main_AppendIncrementalSfM.cpp)
target_link_libraries(openMVG_main_AppendIncrementalSfM
  PRIVATE
    openMVG_system
    openMVG_image
    openMVG_features
    openMVG_sfm
    ${STLPLUS_LIBRARY}
)

set_property(TARGET openMVG_main_AppendIncrementalSfM PROPERTY FOLDER OpenMVG/software)
install(TARGETS openMVG_main_AppendIncrementalSfM DESTINATION bin/)
