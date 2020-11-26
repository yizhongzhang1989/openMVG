//
// Created by v-xinli1 on 11/23/2020.
//

#ifndef OPENMVG_WINDOW_SEQUENTIAL_SFM_HPP
#define OPENMVG_WINDOW_SEQUENTIAL_SFM_HPP

#include <openMVG/sfm/pipelines/sfm_engine.hpp>
#include <openMVG/sfm/sfm.hpp>
#include "openMVG/sfm/pipelines/sequential/sequential_SfM.hpp"

namespace htmlDocument { class htmlDocumentStream; }
namespace { template <typename T> class Histogram; }

namespace openMVG {
namespace sfm {

    struct Features_Provider;
    struct Matches_Provider;

class WindowSequentialSfMReconstructionEngine : public SequentialSfMReconstructionEngine
{
public:

    WindowSequentialSfMReconstructionEngine(
            const SfM_Data & sfm_data,
            const std::string & soutDirectory,
            const std::string & loggingFile = "");

    ~WindowSequentialSfMReconstructionEngine() override;


    virtual bool Process() override;

	bool CreateLocalScene(
        SfM_Data&               local_scene, 
        std::set<IndexT>&       local_scene_viewId, 
        const std::set<IndexT>& add_viewId_cur
    );

    bool BundleAdjustmentWindows(SfM_Data& sfm_data);

    bool badTrackRejectorWindow( SfM_Data& sfm_data, double dPrecision, size_t count = 0);

    double ComputeResidualsHistogram(Histogram<double> * histo);
};

}
}





#endif //OPENMVG_WINDOW_SEQUENTIAL_SFM_HPP
