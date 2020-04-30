// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG_IMU/sfm/pipelines/global/sfm_global_engine_ba.hpp"



#include <iostream>

#ifdef _MSC_VER
#pragma warning( once : 4267 ) //warning C4267: 'argument' : conversion from 'size_t' to 'const int', possible loss of data
#endif

namespace openMVG{
namespace sfm{


GlobalSfMReconstructionEngine_BA::GlobalSfMReconstructionEngine_BA(
  const SfM_Data & sfm_data,
  const std::string & soutDirectory,
  const std::string & sloggingFile)
  : GlobalSfMReconstructionEngine_RelativeMotions(sfm_data,soutDirectory,sloggingFile)
{}

GlobalSfMReconstructionEngine_BA::~GlobalSfMReconstructionEngine_BA()
{
  
}

bool GlobalSfMReconstructionEngine_BA::Process() {

   std::cout<<"/////IMU Global SfM BA/////\n";
  
  if (!Adjust())
  {
    std::cerr << "GlobalSfM:: Non-linear adjustment failure!" << std::endl;
    return false;
  }

  //-- Export statistics about the SfM process
  std::cout << "Structure from Motion statistics.";


  
  std::cout << "-------------------------------" << "\n"
    << "-- View count: " << sfm_data_.GetViews().size() << "\n"
    << "-- Intrinsic count: " << sfm_data_.GetIntrinsics().size() << "\n"
    << "-- Pose count: " << sfm_data_.GetPoses().size() << "\n"
    << "-- Track count: "  << sfm_data_.GetLandmarks().size() << "\n"
    << "-------------------------------" << "\n";
    

  return true;
}

} // namespace sfm
} // namespace openMVG
