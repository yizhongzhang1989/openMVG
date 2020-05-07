// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG_IMU/sfm/pipelines/global/GlobalSfM_priotranslation_averaging.hpp"  //bc

#ifdef OPENMVG_USE_OPENMP
#include <omp.h>
#endif

#include "openMVG/graph/graph.hpp"
#include "openMVG/types.hpp"
#include "openMVG/cameras/Camera_Intrinsics.hpp"
#include "openMVG/cameras/Camera_Pinhole.hpp"
#include "openMVG/linearProgramming/linearProgramming.hpp"
#include "openMVG/multiview/essential.hpp"
#include "openMVG/multiview/translation_averaging_common.hpp"
#include "openMVG_IMU/multiview/translation_averaging_solver.hpp"  //BC
#include "openMVG/robust_estimation/robust_estimator_ACRansac.hpp"
#include "openMVG/sfm/pipelines/global/mutexSet.hpp"
#include "openMVG/sfm/pipelines/global/sfm_global_reindex.hpp"
#include "openMVG/sfm/pipelines/global/triplet_t_ACRansac_kernelAdaptator.hpp"
#include "openMVG/sfm/pipelines/sfm_features_provider.hpp"
#include "openMVG/sfm/pipelines/sfm_matches_provider.hpp"
#include "openMVG/sfm/sfm_data_io.hpp"
#include "openMVG/sfm/sfm_filters.hpp"
#include "openMVG/stl/stl.hpp"
#include "openMVG/system/timer.hpp"

#include <vector>

namespace openMVG{
namespace sfm{

using namespace openMVG::cameras;
using namespace openMVG::geometry;
using namespace openMVG::matching;

/// Use features in normalized camera frames
bool GlobalSfM_PrioTranslation_AveragingSolver::MyRun
(
  ETranslationAveragingMethod eTranslationAveragingMethod,
  openMVG::sfm::SfM_Data & sfm_data,
  const openMVG::sfm::Features_Provider * features_provider,
  const openMVG::sfm::Matches_Provider * matches_provider,
  const Hash_Map<IndexT, Mat3> & map_globalR,
  matching::PairWiseMatches & tripletWise_matches,
  const openMVG::sfm::Matches_Provider * extra_matches_provider
)
{
  // Compute the relative translations and save them to vec_initialRijTijEstimates:
  Compute_translations(
    sfm_data,
    features_provider,
    matches_provider,
    map_globalR,
    tripletWise_matches);
  
  //check whether the extra matches  is filtered completely
  size_t remainextramatches = 0;
  const Pair_Set& pairs = extra_matches_provider->getPairs();
  for (const openMVG::RelativeInfo_Vec & iter : Getrelative_motion())
  {
    for (const relativeInfo & rel : iter)
    {
          if(pairs.count(rel.first)) 
          {
            remainextramatches++;
            break;
          }
    }
  }

  if(!remainextramatches)
  {
    std::cout<<"all extra matches are filtered\n";
    return false;
  }
  std::cout<<"Remain "<<remainextramatches<<" pairs\n";

  const bool b_translation = PrioTranslation_averaging(
    eTranslationAveragingMethod,
    sfm_data,
    map_globalR);

  // Filter matches to keep only them link to view that have valid poses
  // (necessary since multiple components exists before translation averaging)
  std::set<IndexT> valid_view_ids;
  for (const auto & view : sfm_data.GetViews())
  {
    if (sfm_data.IsPoseAndIntrinsicDefined(view.second.get()))
      valid_view_ids.insert(view.first);
  }
  KeepOnlyReferencedElement(valid_view_ids, tripletWise_matches);

  return b_translation;
}

bool GlobalSfM_PrioTranslation_AveragingSolver::PrioTranslation_averaging(
  ETranslationAveragingMethod eTranslationAveragingMethod,
  sfm::SfM_Data & sfm_data,
  const Hash_Map<IndexT, Mat3> & map_globalR)
{
  //-------------------
  //-- GLOBAL TRANSLATIONS ESTIMATION from initial triplets t_ij guess
  //-------------------

  // Keep the largest Biedge connected component graph of relative translations
  Pair_Set pairs;
  for (const openMVG::RelativeInfo_Vec & iter : vec_relative_motion_)
  {
    for (const relativeInfo & rel : iter)
    {
      pairs.insert(rel.first);
    }
  }
  const std::set<IndexT> set_remainingIds =
    openMVG::graph::CleanGraph_KeepLargestBiEdge_Nodes<Pair_Set, IndexT>(pairs);
  KeepOnlyReferencedElement(set_remainingIds, vec_relative_motion_);

  {
    const std::set<IndexT> index = getIndexT(vec_relative_motion_);

    const size_t iNview = index.size();
    std::cout << "\n-------------------------------" << "\n"
      << " Global translations computation: " << "\n"
      << "   - Ready to compute " << iNview << " global translations." << "\n"
      << "     from #relative translations: " << vec_relative_motion_.size()*3 << std::endl;

    if (iNview < 3)
    {
      // Too tiny image set to perform motion averaging
      return false;
    }
    //-- Update initial estimates from [minId,maxId] to range [0->Ncam]
    std::vector<RelativeInfo_Vec> vec_relative_motion_cpy = vec_relative_motion_;
    const Pair_Set pairs = getPairs(vec_relative_motion_cpy);
    Hash_Map<IndexT,IndexT> reindex_forward, reindex_backward;
    reindex(pairs, reindex_forward, reindex_backward);
    for (openMVG::RelativeInfo_Vec & iter : vec_relative_motion_cpy)
    {
      for (relativeInfo & it : iter)
      {
        it.first = Pair(reindex_forward[it.first.first], reindex_forward[it.first.second]);
      }
    }

    openMVG::system::Timer timerLP_translation;

    switch (eTranslationAveragingMethod)
    {
      case TRANSLATION_AVERAGING_L1:
      {
        double gamma = -1.0;
        std::vector<double> vec_solution;
        {
          vec_solution.resize(iNview*3 + vec_relative_motion_cpy.size() + 1);
          using namespace openMVG::linearProgramming;
          OSI_CLP_SolverWrapper solverLP(vec_solution.size());

          lInfinityCV::Tifromtij_ConstraintBuilder cstBuilder(vec_relative_motion_cpy);

          LP_Constraints_Sparse constraint;
          //-- Setup constraint and solver
          cstBuilder.Build(constraint);
          solverLP.setup(constraint);
          //--
          // Solving
          const bool bFeasible = solverLP.solve();
          std::cout << " \n Feasibility " << bFeasible << std::endl;
          //--
          if (bFeasible)  {
            solverLP.getSolution(vec_solution);
            gamma = vec_solution[vec_solution.size()-1];
          }
          else  {
            std::cerr << "Compute global translations: failed" << std::endl;
            return false;
          }
        }

        const double timeLP_translation = timerLP_translation.elapsed();
        //-- Export triplet statistics:
        {

          std::ostringstream os;
          os << "Translation fusion statistics.";
          os.str("");
          os << "-------------------------------" << "\n"
            << "-- #relative estimates: " << vec_relative_motion_cpy.size()
            << " converge with gamma: " << gamma << ".\n"
            << " timing (s): " << timeLP_translation << ".\n"
            << "-------------------------------" << "\n";
          std::cout << os.str() << std::endl;
        }

        std::cout << "Found solution:\n";
        std::copy(vec_solution.begin(), vec_solution.end(), std::ostream_iterator<double>(std::cout, " "));

        std::vector<double> vec_camTranslation(iNview*3,0);
        std::copy(&vec_solution[0], &vec_solution[iNview*3], &vec_camTranslation[0]);

        std::vector<double> vec_camRelLambdas(&vec_solution[iNview*3], &vec_solution[iNview*3 + vec_relative_motion_cpy.size()]);
        std::cout << "\ncam position: " << std::endl;
        std::copy(vec_camTranslation.begin(), vec_camTranslation.end(), std::ostream_iterator<double>(std::cout, " "));
        std::cout << "\ncam Lambdas: " << std::endl;
        std::copy(vec_camRelLambdas.begin(), vec_camRelLambdas.end(), std::ostream_iterator<double>(std::cout, " "));
        std::cout << std::endl;

        // Update the view poses according the found camera centers
        for (size_t i = 0; i < iNview; ++i)
        {
          const Vec3 t(vec_camTranslation[i*3], vec_camTranslation[i*3+1], vec_camTranslation[i*3+2]);
          const IndexT pose_id = reindex_backward[i];
          const Mat3 & Ri = map_globalR.at(pose_id);
          sfm_data.poses[pose_id] = Pose3(Ri, -Ri.transpose()*t);
        }
      }
      break;

      case TRANSLATION_AVERAGING_SOFTL1:
      {
        std::vector<Vec3> vec_translations;
		///////BC START////////
        std::vector<Vec3> prio_translations;
		prio_translations.reserve(sfm_data.GetPoses().size());
        
        for(const auto& pose_item:sfm_data.GetPoses())
        {
            prio_translations.emplace_back(pose_item.second.translation());
        }
        
        if (!solve_priotranslations_problem_softl1(
          vec_relative_motion_cpy, vec_translations,prio_translations))
        {
          std::cerr << "Compute global translations: failed" << std::endl;
          return false;
        }
		///////BC END//////////
        // A valid solution was found:
        // - Update the view poses according the found camera translations
        for (size_t i = 0; i < iNview; ++i)
        {
          const Vec3 & t = vec_translations[i];
          const IndexT pose_id = reindex_backward[i];
          const Mat3 & Ri = map_globalR.at(pose_id);
          sfm_data.poses[pose_id] = Pose3(Ri, -Ri.transpose()*t);
        }
      }
      break;

      case TRANSLATION_AVERAGING_L2_DISTANCE_CHORDAL:
      {
        std::vector<int> vec_edges;
        vec_edges.reserve(vec_relative_motion_cpy.size() * 2);
        std::vector<double> vec_poses;
        vec_poses.reserve(vec_relative_motion_cpy.size() * 3);
        std::vector<double> vec_weights;
        vec_weights.reserve(vec_relative_motion_cpy.size());

        for (const openMVG::RelativeInfo_Vec & iter : vec_relative_motion_cpy)
        {
          for (const relativeInfo & rel : iter)
          {
            vec_edges.push_back(rel.first.first);
            vec_edges.push_back(rel.first.second);
            // Since index have been remapped
            // (use the backward indexing to retrieve the second global rotation)
            const IndexT secondId = reindex_backward[rel.first.second];
            const View * view = sfm_data.views.at(secondId).get();
            const Mat3 & Ri = map_globalR.at(view->id_pose);
            const Vec3 direction = -(Ri.transpose() * rel.second.second.normalized());

            vec_poses.push_back(direction(0));
            vec_poses.push_back(direction(1));
            vec_poses.push_back(direction(2));

            vec_weights.push_back(1.0);
          }
        }

        const double function_tolerance = 1e-7, parameter_tolerance = 1e-8;
        const int max_iterations = 500;

        const double loss_width = 0.0; // No loss in order to compare with TRANSLATION_AVERAGING_L1

        std::vector<double> X(iNview*3, 0.0);
        if (!solve_translations_problem_l2_chordal(
          &vec_edges[0],
          &vec_poses[0],
          &vec_weights[0],
          vec_relative_motion_cpy.size()*3,
          loss_width,
          &X[0],
          function_tolerance,
          parameter_tolerance,
          max_iterations))  {
            std::cerr << "Compute global translations: failed" << std::endl;
            return false;
        }

        // Update camera center for each view
        for (size_t i = 0; i < iNview; ++i)
        {
          const Vec3 C(X[i*3], X[i*3+1], X[i*3+2]);
          const IndexT pose_id = reindex_backward[i]; // undo the reindexing
          const Mat3 & Ri = map_globalR.at(pose_id);
          sfm_data.poses[pose_id] = Pose3(Ri, C);
        }
      }
      break;
      default:
      {
        std::cerr << "Unknown translation averaging method" << std::endl;
        return false;
      }
    }
  }
  return true;
}



} // namespace sfm
} // namespace openMVG
