
add_executable(openMVG_main_IncrementalSfM_window
        main_IncrementalSfM_window.cpp
        window_sequential_SfM.cpp
        )
target_link_libraries(openMVG_main_IncrementalSfM_window
        PRIVATE
        openMVG_system
        openMVG_image
        openMVG_features
        openMVG_sfm
        ${STLPLUS_LIBRARY}
        )

target_include_directories(openMVG_main_IncrementalSfM_window PRIVATE ${CERES_INCLUDE_DIRS})
set_property(TARGET openMVG_main_IncrementalSfM_window PROPERTY FOLDER OpenMVG/software/Window_SfM)
install(TARGETS openMVG_main_IncrementalSfM_window DESTINATION bin/)