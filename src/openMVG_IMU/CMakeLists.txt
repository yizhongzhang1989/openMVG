add_subdirectory(sfm)
add_subdirectory(multiview)
add_subdirectory(matching_image_collection)
add_subdirectory(matching)
add_subdirectory(utils)
#INSTALL RULES
install(
  DIRECTORY .
  DESTINATION include/openMVG
  COMPONENT headers
  FILES_MATCHING PATTERN "*.hpp" PATTERN "*.h"
)
