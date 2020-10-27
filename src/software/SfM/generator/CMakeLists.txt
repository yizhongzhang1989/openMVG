add_library(simulation_generator SimulationGenerator.cpp)

add_executable(openMVG_main_GenerateSimulationData main_GenerateSimulationData.cpp)
target_link_libraries(openMVG_main_GenerateSimulationData
    PRIVATE
        openMVG_system
        openMVG_image
        openMVG_features
        openMVG_sfm
        ${STLPLUS_LIBRARY}
        vlsift
        simulation_generator)

set_property(TARGET openMVG_main_GenerateSimulationData PROPERTY FOLDER OpenMVG/software)
install(TARGETS openMVG_main_GenerateSimulationData DESTINATION bin/)