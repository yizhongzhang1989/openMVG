add_executable(openMVG_main_simu_gen main_simu_gen.cpp)
target_link_libraries(openMVG_main_simu_gen
        PRIVATE
        openMVG_system
        openMVG_image
        openMVG_features
        openMVG_sfm
        ${STLPLUS_LIBRARY}
        vlsift)

set_property(TARGET openMVG_main_simu_gen PROPERTY FOLDER OpenMVG/software)
install(TARGETS openMVG_main_simu_gen DESTINATION bin/)

add_executable(main_simu_gen_time main_simu_gen_time.cpp)