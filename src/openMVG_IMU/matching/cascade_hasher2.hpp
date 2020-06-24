// This file is part of OpenMVG_IMU , a branch of OpenMVG
// Author: Bao Chong
// Date:2020/06

#ifndef OPENMVG_IMU_MATCHING_CASCADE_2_HASHER_HPP
#define OPENMVG_IMU_MATCHING_CASCADE_2_HASHER_HPP

//------------------
//-- Bibliography --
//------------------
//- [1] "Fast and Accurate Image Matching with Cascade Hashing for 3D Reconstruction"
//- Authors: Jian Cheng, Cong Leng, Jiaxiang Wu, Hainan Cui, Hanqing Lu.
//- Date: 2014.
//- Conference: CVPR.
//
// This implementation is based on the Theia library implementation.
//
// Update compare to the initial paper [1] and initial author code:
// - hashing projection is made by using Eigen to use vectorization (Theia)
// - replace the BoxMuller random number generation by C++ 11 random number generation (OpenMVG)
// - this implementation can support various descriptor length and internal type (OpenMVG)
// -  SIFT, SURF, ... all scalar based descriptor
//

// Copyright (C) 2014 The Regents of the University of California (Regents).
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
//
//     * Neither the name of The Regents or University of California nor the
//       names of its contributors may be used to endorse or promote products
//       derived from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDERS OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//
// Please contact the author of this library if you have any questions.
// Author: Chris Sweeney (cmsweeney@cs.ucsb.edu)


#include <cmath>
#include <random>
#include <utility>
#include <vector>

#include "openMVG/matching/indMatch.hpp"
#include "openMVG_IMU/matching/cascade_hasher.hpp"

namespace openMVG {
namespace matching {



// The class inherit the `CascadeHasher_General` in order to 
// provide a matching interface that can only perform sift matching 
// on user specified features of first view.
class CascadeHasher2:public CascadeHasher_General {

public:
  CascadeHasher2() = default;
  template <typename MatrixT, typename DistanceType>
  void Match_HashedDescriptions_Specifiedfeat1
  (
	  const HashedDescriptions& hashed_descriptions1,
	  const MatrixT & descriptions1,
	  const HashedDescriptions& hashed_descriptions2,
	  const MatrixT & descriptions2,
	  IndMatches * pvec_indices,
	  std::vector<DistanceType> * pvec_distances,
	  const std::set<IndexT>& specified_feat1_ids,  //user specified feature ids
	  const int NN = 2
  )
  {
	  {
		  using MetricT = L2<typename MatrixT::Scalar>;
		  MetricT metric;

		  static const int kNumTopCandidates = 10;

		  // Preallocate the candidate descriptors container.
		  std::vector<int> candidate_descriptors;
		  candidate_descriptors.reserve(hashed_descriptions2.hashed_desc.size());

		  // Preallocated hamming distances. Each column indicates the hamming distance
		  // and the rows collect the descriptor ids with that
		  // distance. num_descriptors_with_hamming_distance keeps track of how many
		  // descriptors have that distance.
		  Eigen::MatrixXi candidate_hamming_distances(
			  hashed_descriptions2.hashed_desc.size(), nb_hash_code_ + 1);
		  Eigen::VectorXi num_descriptors_with_hamming_distance(nb_hash_code_ + 1);

		  // Preallocate the container for keeping euclidean distances.
		  std::vector<std::pair<DistanceType, int>> candidate_euclidean_distances;
		  candidate_euclidean_distances.reserve(kNumTopCandidates);

		  // A preallocated vector to determine if we have already used a particular
		  // feature for matching (i.e., prevents duplicates).
		  std::vector<bool> used_descriptor(hashed_descriptions2.hashed_desc.size());

		  using HammingMetricType = matching::Hamming<stl::dynamic_bitset::BlockType>;
		  static const HammingMetricType metricH = {};
		  for (int i = 0; i < hashed_descriptions1.hashed_desc.size(); ++i)
		  {
			  if (!specified_feat1_ids.count(i)) continue;   //bc
			  candidate_descriptors.clear();
			  num_descriptors_with_hamming_distance.setZero();
			  candidate_euclidean_distances.clear();

			  const auto& hashed_desc = hashed_descriptions1.hashed_desc[i];

			  // Accumulate all descriptors in each bucket group that are in the same
			  // bucket id as the query descriptor.
			  for (int j = 0; j < nb_bucket_groups_; ++j)
			  {
				  const uint16_t bucket_id = hashed_desc.bucket_ids[j];
				  for (const auto& feature_id : hashed_descriptions2.buckets[j][bucket_id])
				  {
					  candidate_descriptors.emplace_back(feature_id);
					  used_descriptor[feature_id] = false;
				  }
			  }

			  // Skip matching this descriptor if there are not at least NN candidates.
			  if (candidate_descriptors.size() <= NN)
			  {
				  continue;
			  }

			  // Compute the hamming distance of all candidates based on the comp hash
			  // code. Put the descriptors into buckets corresponding to their hamming
			  // distance.
			  for (const int candidate_id : candidate_descriptors)
			  {
				  if (!used_descriptor[candidate_id]) // avoid selecting the same candidate multiple times
				  {
					  used_descriptor[candidate_id] = true;

					  const HammingMetricType::ResultType hamming_distance = metricH(
						  hashed_desc.hash_code.data(),
						  hashed_descriptions2.hashed_desc[candidate_id].hash_code.data(),
						  hashed_desc.hash_code.num_blocks());
					  candidate_hamming_distances(
						  num_descriptors_with_hamming_distance(hamming_distance)++,
						  hamming_distance) = candidate_id;
				  }
			  }

			  // Compute the euclidean distance of the k descriptors with the best hamming
			  // distance.
			  candidate_euclidean_distances.reserve(kNumTopCandidates);
			  for (int j = 0; j < candidate_hamming_distances.cols() &&
				  (candidate_euclidean_distances.size() < kNumTopCandidates); ++j)
			  {
				  for (int k = 0; k < num_descriptors_with_hamming_distance(j) &&
					  (candidate_euclidean_distances.size() < kNumTopCandidates); ++k)
				  {
					  const int candidate_id = candidate_hamming_distances(k, j);
					  const DistanceType distance = metric(
						  descriptions2.row(candidate_id).data(),
						  descriptions1.row(i).data(),
						  descriptions1.cols());

					  candidate_euclidean_distances.emplace_back(distance, candidate_id);
				  }
			  }

			  // Assert that each query is having at least NN retrieved neighbors
			  if (candidate_euclidean_distances.size() >= NN)
			  {
				  // Find the top NN candidates based on euclidean distance.
				  std::partial_sort(candidate_euclidean_distances.begin(),
					  candidate_euclidean_distances.begin() + NN,
					  candidate_euclidean_distances.end());
				  // save resulting neighbors
				  for (int l = 0; l < NN; ++l)
				  {
					  pvec_distances->emplace_back(candidate_euclidean_distances[l].first);
					  pvec_indices->emplace_back(IndMatch(i, candidate_euclidean_distances[l].second));
				  }
			  }
			  //else -> too few candidates... (save no one)
		  }
	  }
  }
  

};



}  // namespace matching
}  // namespace openMVG

#endif // OPENMVG_MATCHING_CASCADE_HASHER_HPP
