find_package(Pangolin REQUIRED)

include_directories(${EIGEN_INCLUDE_DIRS} ${Pangolin_INCLUDE_DIRS})

add_library(simulation_generator SimulationGenerator.cpp)
set_property(TARGET simulation_generator PROPERTY FOLDER OpenMVG/software)
install(TARGETS simulation_generator DESTINATION lib/)

add_executable(openMVG_main_GenerateSimulationData main_GenerateSimulationData.cpp)
target_link_libraries(openMVG_main_GenerateSimulationData
    PRIVATE
        openMVG_system
        openMVG_image
        openMVG_features
        openMVG_sfm
        ${STLPLUS_LIBRARY}
        vlsift
        simulation_generator
        ${Pangolin_LIBRARIES})

set_property(TARGET openMVG_main_GenerateSimulationData PROPERTY FOLDER OpenMVG/software)
install(TARGETS openMVG_main_GenerateSimulationData DESTINATION bin/)
