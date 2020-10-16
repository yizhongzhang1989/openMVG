add_executable(openMVG_main_VISfM main_VISfM.cpp
        VISfM_src/sequential_VISfM.cpp VISfM_src/sequential_VISfM.hpp VISfM_src/visfm_data_BA_ceres.cpp VISfM_src/visfm_data_BA_ceres.hpp VISfM_src/VISfM_ceres_facotr.cpp VISfM_src/VISfM_ceres_facotr.h VISfM_src/VISfM_ceres_param.cpp VISfM_src/VISfM_ceres_param.h)
target_link_libraries(openMVG_main_VISfM
        PRIVATE
        openMVG_system
        openMVG_image
        openMVG_features
        openMVG_sfm
        ${STLPLUS_LIBRARY}
        )

set_property(TARGET openMVG_main_VISfM PROPERTY FOLDER OpenMVG/software/VISfM)
install(TARGETS openMVG_main_VISfM DESTINATION bin/)