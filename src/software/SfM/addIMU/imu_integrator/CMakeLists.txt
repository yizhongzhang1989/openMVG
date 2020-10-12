#cmake_minimum_required(VERSION 2.8)
#project(imu_integrator)

#find_package(Eigen3 REQUIRED)

include_directories(
    ${EIGEN3_INCLUDE_DIRS}
#    ${PROJECT_SOURCE_DIR}
)

#set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/lib)

aux_source_directory(./imu_integrator/ IMU_INTE_SRCS)
add_library(imu_integrator SHARED ${IMU_INTE_SRCS})

set_property(TARGET imu_integrator PROPERTY FOLDER OpenMVG/software)
install(TARGETS imu_integrator DESTINATION lib/)

#add_executable(test ./test.cpp)
#target_link_libraries(test ${PROJECT_NAME})

