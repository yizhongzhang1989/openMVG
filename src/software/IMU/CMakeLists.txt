
###
# Intrinsic image analysis and SfM_Data container initialization
###

# Installation rules


###
# Add executable that computes:
# - openMVG_main_ComputeFeatures: features and descriptors
# - openMVG_main_ComputeMatches: putative matches + geometric filtered matches
###


# Installation rules


###
# SfM Pipelines
# - Incremental/Sequential
# - Global
###
# - Convert SfM_Data format (from one to another)
# - Pair filtering (from frustum)
# - Compute structure from known camera poses
# - Compute structure color
###


add_executable(openMVG_IMU_main_GlobalSfM_motion_averaging main_GlobalSfM_motion_averaging.cpp)
target_link_libraries(openMVG_IMU_main_GlobalSfM_motion_averaging
  PRIVATE
    openMVG_IMU_sfm
    openMVG_system
    openMVG_image
    openMVG_features
    openMVG_sfm
    ${STLPLUS_LIBRARY}
)

# Installation rules
set_property(TARGET openMVG_IMU_main_GlobalSfM_motion_averaging PROPERTY FOLDER OpenMVG/software/imu)
install(TARGETS openMVG_IMU_main_GlobalSfM_motion_averaging DESTINATION bin/)

#Bundle Adjustment in Global SfM
add_executable(openMVG_IMU_main_GlobalSfM_ba main_GlobalSfM_ba.cpp)
target_link_libraries(openMVG_IMU_main_GlobalSfM_ba
  PRIVATE
    openMVG_IMU_sfm
    openMVG_system
    openMVG_sfm
    ${STLPLUS_LIBRARY}
)

# Installation rules
set_property(TARGET openMVG_IMU_main_GlobalSfM_ba PROPERTY FOLDER OpenMVG/software/imu)
install(TARGETS openMVG_IMU_main_GlobalSfM_ba DESTINATION bin/)

#Pose Augmentation in Global SfM
add_executable(openMVG_IMU_main_GlobalSfM_pose_augmentation main_GlobalSfM_pose_augmentation.cpp)
target_link_libraries(openMVG_IMU_main_GlobalSfM_pose_augmentation
  PRIVATE
    openMVG_IMU_sfm
    openMVG_system
    openMVG_sfm
    ${STLPLUS_LIBRARY}
)

# Installation rules
set_property(TARGET openMVG_IMU_main_GlobalSfM_pose_augmentation PROPERTY FOLDER OpenMVG/software/imu)
install(TARGETS openMVG_IMU_main_GlobalSfM_pose_augmentation DESTINATION bin/)

#Iterative Pose Augmentation in Global SfM
add_executable(openMVG_IMU_main_GlobalSfM_iterative_pose_augmentation main_GlobalSfM_iterative_pose_augmentation.cpp)
target_link_libraries(openMVG_IMU_main_GlobalSfM_iterative_pose_augmentation
  PRIVATE
    openMVG_IMU_sfm
    openMVG_system
    openMVG_sfm
    ${STLPLUS_LIBRARY}
)

# Installation rules
set_property(TARGET openMVG_IMU_main_GlobalSfM_iterative_pose_augmentation PROPERTY FOLDER OpenMVG/software/imu)
install(TARGETS openMVG_IMU_main_GlobalSfM_iterative_pose_augmentation DESTINATION bin/)

###
# SfM Align tools
###
add_executable(openMVG_IMU_main_TrajectoryAlign main_TrajectoryAlign.cpp)
target_link_libraries(openMVG_IMU_main_TrajectoryAlign
  openMVG_features
  openMVG_image
  openMVG_sfm
  openMVG_system
  ${STLPLUS_LIBRARY}
)

#installation rules
set_property(TARGET openMVG_IMU_main_TrajectoryAlign PROPERTY FOLDER OpenMVG/software/imu)
install(TARGETS openMVG_IMU_main_TrajectoryAlign DESTINATION bin/)
###
# Iterative sequential sfm
###
add_executable(openMVG_IMU_main_Iterative_IncrementalSfM main_Iterative_IncrementalSfM.cpp)
target_link_libraries(openMVG_IMU_main_Iterative_IncrementalSfM
  PRIVATE
    openMVG_system
    openMVG_image
    openMVG_features
    openMVG_sfm
    ${STLPLUS_LIBRARY}
)

set_property(TARGET openMVG_IMU_main_Iterative_IncrementalSfM PROPERTY FOLDER OpenMVG/software/imu)
install(TARGETS openMVG_IMU_main_Iterative_IncrementalSfM DESTINATION bin/)


