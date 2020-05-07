add_subdirectory(sfm)
add_subdirectory(multiview)

#INSTALL RULES
install(
  DIRECTORY .
  DESTINATION include/openMVG
  COMPONENT headers
  FILES_MATCHING PATTERN "*.hpp" PATTERN "*.h"
)
