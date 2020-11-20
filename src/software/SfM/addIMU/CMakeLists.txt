add_subdirectory(imu_integrator)

# include_directories(${CERES_INCLUDE_DIRS})

# add_executable(openMVG_main_IncrementalSfM_IMU 
#     main_IncrementalSfM_IMU.cpp 
#     sequential_SfM_IMU.cpp
#     sfm_imu.cpp)

# target_link_libraries(openMVG_main_IncrementalSfM_IMU
#   PRIVATE
#     openMVG_system
#     openMVG_image
#     openMVG_features
#     openMVG_sfm
#     ${STLPLUS_LIBRARY}
#     imu_integrator
# )

# add_executable(openMVG_main_SfMInit_ImageListingFromVideo main_SfMInit_ImageListingFromVideo.cpp)
# target_link_libraries(openMVG_main_SfMInit_ImageListingFromVideo
#   PRIVATE
#     openMVG_system
#     openMVG_image
#     openMVG_features
#     openMVG_sfm
#     openMVG_exif
#     ${STLPLUS_LIBRARY}
# )

# set_property(TARGET openMVG_main_IncrementalSfM_IMU PROPERTY FOLDER OpenMVG/software)
# install(TARGETS openMVG_main_IncrementalSfM_IMU DESTINATION bin/)
# set_property(TARGET openMVG_main_SfMInit_ImageListingFromVideo PROPERTY FOLDER OpenMVG/software)
# install(TARGETS openMVG_main_SfMInit_ImageListingFromVideo DESTINATION bin/)
