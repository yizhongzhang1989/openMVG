
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


add_executable(openMVG_IMU_main_GlobalSfM main_GlobalSfM.cpp)
target_link_libraries(openMVG_IMU_main_GlobalSfM
  PRIVATE
    openMVG_IMU_sfm
    openMVG_system
    openMVG_image
    openMVG_features
    openMVG_sfm
    ${STLPLUS_LIBRARY}
)

# Installation rules
set_property(TARGET openMVG_IMU_main_GlobalSfM PROPERTY FOLDER OpenMVG/software/imu)
install(TARGETS openMVG_IMU_main_GlobalSfM DESTINATION bin/)

# SplitMatchFileIntoMatchFiles


# Installation rules


# -
# Use openCV Feature detector to describe the images
# -


###
# SfM Research tools
#
# - Quality comparison against a GT camera path (MultiView Evaluation dataset)
###


#installation rules


##
# Export tutorial with valid path thanks to a configuration file
##



