// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2015 Pierre Moulon.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_SFM_SFM_DATA_IMU_BA_CERES_HPP
#define OPENMVG_SFM_SFM_DATA_IMU_BA_CERES_HPP


#include <openMVG/sfm/sfm_data_BA_ceres.hpp>

namespace ceres { class CostFunction; }
namespace openMVG { namespace cameras { struct IntrinsicBase; } }
namespace openMVG { namespace sfm { struct SfM_Data; } }

namespace openMVG {
namespace sfm {

class Bundle_Adjustment_IMU_Ceres : public Bundle_Adjustment_Ceres
{
public:
    explicit Bundle_Adjustment_IMU_Ceres
        (
            const Bundle_Adjustment_Ceres::BA_Ceres_options & options =
            std::move(BA_Ceres_options())
        ):Bundle_Adjustment_Ceres(options){}
//    bool Adjust
//        (
//            // the SfM scene to refine
//            sfm::SfM_Data & sfm_data,
//            // tell which parameter needs to be adjusted
//            const Optimize_Options & options
//        ) override;
};

} // namespace sfm
} // namespace openMVG

#endif // OPENMVG_SFM_SFM_DATA_IMU_BA_CERES_HPP
