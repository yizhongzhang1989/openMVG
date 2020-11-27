add_subdirectory(VISfM_src)

add_executable(openMVG_main_VISfM main_VISfM.cpp)
#        VISfM_src/sequential_VISfM.cpp VISfM_src/sequential_VISfM.hpp VISfM_src/visfm_data_BA_ceres.cpp VISfM_src/visfm_data_BA_ceres.hpp VISfM_src/VISfM_ceres_facotr.cpp VISfM_src/VISfM_ceres_facotr.hpp VISfM_src/VISfM_ceres_param.cpp VISfM_src/VISfM_ceres_param.hpp VISfM_src/Utility.hpp VISfM_src/VI_static_Parm.hpp VISfM_src/IMU_InteBase.hpp VISfM_src/VI_static_Parm.cpp)
target_link_libraries(openMVG_main_VISfM
        PRIVATE
        openMVG_vi_sfm
        openMVG_system
        openMVG_image
        openMVG_features
        openMVG_sfm
        ${STLPLUS_LIBRARY}
        )
target_include_directories(openMVG_main_VISfM PRIVATE ${CERES_INCLUDE_DIRS})
set_property(TARGET openMVG_main_VISfM PROPERTY FOLDER OpenMVG/software/VISfM)
install(TARGETS openMVG_main_VISfM DESTINATION bin/)

add_executable(openMVG_main_OnlyScaleEstiamtion OnlyScaleEstiamtion.cpp)
target_link_libraries(openMVG_main_OnlyScaleEstiamtion
        PRIVATE
        openMVG_vi_sfm
        openMVG_system
        openMVG_image
        openMVG_features
        openMVG_sfm
        ${STLPLUS_LIBRARY}
        )
target_include_directories(openMVG_main_OnlyScaleEstiamtion PRIVATE ${CERES_INCLUDE_DIRS})
set_property(TARGET openMVG_main_OnlyScaleEstiamtion PROPERTY FOLDER OpenMVG/software/VISfM)
install(TARGETS openMVG_main_OnlyScaleEstiamtion DESTINATION bin/)


add_executable(openMVG_VISfM_test VISfM_test.cpp)
#        VISfM_src/sequential_VISfM.cpp VISfM_src/sequential_VISfM.hpp VISfM_src/visfm_data_BA_ceres.cpp VISfM_src/visfm_data_BA_ceres.hpp VISfM_src/VISfM_ceres_facotr.cpp VISfM_src/VISfM_ceres_facotr.hpp VISfM_src/VISfM_ceres_param.cpp VISfM_src/VISfM_ceres_param.hpp VISfM_src/Utility.hpp VISfM_src/VI_static_Parm.hpp VISfM_src/IMU_InteBase.hpp VISfM_src/VI_static_Parm.cpp)
target_link_libraries(openMVG_VISfM_test
        PRIVATE
        openMVG_vi_sfm
        openMVG_system
        openMVG_image
        openMVG_features
        openMVG_sfm
        ${STLPLUS_LIBRARY}
        )
target_include_directories(openMVG_VISfM_test PRIVATE ${CERES_INCLUDE_DIRS})
set_property(TARGET openMVG_VISfM_test PROPERTY FOLDER OpenMVG/software/VISfM)
install(TARGETS openMVG_VISfM_test DESTINATION bin/)

add_executable(main_window_VISfM main_window_VISfM.cpp)
#        VISfM_src/sequential_VISfM.cpp VISfM_src/sequential_VISfM.hpp VISfM_src/visfm_data_BA_ceres.cpp VISfM_src/visfm_data_BA_ceres.hpp VISfM_src/VISfM_ceres_facotr.cpp VISfM_src/VISfM_ceres_facotr.hpp VISfM_src/VISfM_ceres_param.cpp VISfM_src/VISfM_ceres_param.hpp VISfM_src/Utility.hpp VISfM_src/VI_static_Parm.hpp VISfM_src/IMU_InteBase.hpp VISfM_src/VI_static_Parm.cpp)
target_link_libraries(main_window_VISfM
        PRIVATE
        openMVG_vi_sfm
        openMVG_system
        openMVG_image
        openMVG_features
        openMVG_sfm
        ${STLPLUS_LIBRARY}
        )
target_include_directories(main_window_VISfM PRIVATE ${CERES_INCLUDE_DIRS})
set_property(TARGET main_window_VISfM PROPERTY FOLDER OpenMVG/software/VISfM)
install(TARGETS main_window_VISfM DESTINATION bin/)


add_executable(test_imu_factor test_imu_factor.cpp)
#        VISfM_src/sequential_VISfM.cpp VISfM_src/sequential_VISfM.hpp VISfM_src/visfm_data_BA_ceres.cpp VISfM_src/visfm_data_BA_ceres.hpp VISfM_src/VISfM_ceres_facotr.cpp VISfM_src/VISfM_ceres_facotr.hpp VISfM_src/VISfM_ceres_param.cpp VISfM_src/VISfM_ceres_param.hpp VISfM_src/Utility.hpp VISfM_src/VI_static_Parm.hpp VISfM_src/IMU_InteBase.hpp VISfM_src/VI_static_Parm.cpp)
target_link_libraries(test_imu_factor
        PRIVATE
        openMVG_vi_sfm
        openMVG_system
        openMVG_image
        openMVG_features
        openMVG_sfm
        ${STLPLUS_LIBRARY}
        )
target_include_directories(test_imu_factor PRIVATE ${CERES_INCLUDE_DIRS})
set_property(TARGET test_imu_factor PROPERTY FOLDER OpenMVG/software/VISfM)
install(TARGETS test_imu_factor DESTINATION bin/)


add_executable(openMVG_main_save_pose_point save_pose_point.cpp)
target_link_libraries(openMVG_main_save_pose_point
        PRIVATE
        openMVG_system
        openMVG_image
        openMVG_features
        openMVG_sfm
        ${STLPLUS_LIBRARY}
        )
set_property(TARGET openMVG_main_save_pose_point PROPERTY FOLDER OpenMVG/software/VISfM)
install(TARGETS openMVG_main_save_pose_point DESTINATION bin/)

add_executable(openMVG_main_save_tum save_tum.cpp)
target_link_libraries(openMVG_main_save_tum
        PRIVATE
        openMVG_system
        openMVG_image
        openMVG_features
        openMVG_sfm
        ${STLPLUS_LIBRARY}
        )

add_executable(change_time_xinli change_time_xinli.cpp)

add_executable(GetSimuTime GetSimuTime.cpp)

set_property(TARGET openMVG_main_save_tum PROPERTY FOLDER OpenMVG/software/VISfM)
install(TARGETS openMVG_main_save_tum DESTINATION bin/)