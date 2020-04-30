// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_SFM_GLOBAL_ENGINE_BA_HPP
#define OPENMVG_SFM_GLOBAL_ENGINE_BA_HPP

#include <memory>
#include <string>

#include "openMVG/sfm/pipelines/sfm_engine.hpp"
#include "openMVG/sfm/pipelines/global/sfm_global_engine_relative_motions.hpp"

namespace htmlDocument { class htmlDocumentStream; }
namespace openMVG { namespace matching { struct PairWiseMatches; } }
namespace openMVG { namespace sfm { struct Features_Provider; } }
namespace openMVG { namespace sfm { struct Matches_Provider; } }
namespace openMVG { namespace sfm { struct SfM_Data; } }

namespace openMVG{
namespace sfm{

/// Global SfM Pipeline Reconstruction Engine.
/// - Method: Global Fusion of Relative Motions.
class GlobalSfMReconstructionEngine_BA : public GlobalSfMReconstructionEngine_RelativeMotions
{
public:

  GlobalSfMReconstructionEngine_BA(
    const SfM_Data & sfm_data,
    const std::string & soutDirectory,
    const std::string & loggingFile = "");

  ~GlobalSfMReconstructionEngine_BA();

  bool Process();

};

} // namespace sfm
} // namespace openMVG

#endif // OPENMVG_SFM_GLOBAL_ENGINE_RELATIVE_MOTIONS_HPP
