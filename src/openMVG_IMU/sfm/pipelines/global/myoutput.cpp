// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG_IMU/sfm/pipelines/global/myoutput.hpp"

#include <iostream>
#include <fstream>
#include "openMVG/stl/stl.hpp"



#ifdef _MSC_VER
#pragma warning( once : 4267 ) //warning C4267: 'argument' : conversion from 'size_t' to 'const int', possible loss of data
#endif

namespace openMVG{
namespace sfm{



bool Output_trajectory(std::string filename,const SfM_Data& sfm_data_) {

	//save the raw pose
  std::cout<<"Output_trajectory\n";
	std::ofstream trajectory_file(filename);
	trajectory_file << "viewid,R00,R01,R02,R10,R11,R12,R20,R21,R22,C0,C1,C2\n";
	for (const auto& view_item : sfm_data_.GetViews())
	{
		if (!sfm_data_.GetPoses().count(view_item.first))
		{
			std::cout << "View " << view_item.first << " is not calibrated\n";
			continue;
		}
		const geometry::Pose3& pose = sfm_data_.GetPoseOrDie(view_item.second.get());
		trajectory_file << view_item.first << ",";
		const Mat3& rotation = pose.rotation();
		const Vec3& center = pose.center();
		trajectory_file << rotation(0, 0) << ","
			<< rotation(0, 1) << ","
			<< rotation(0, 2) << ","
			<< rotation(1, 0) << ","
			<< rotation(1, 1) << ","
			<< rotation(1, 2) << ","
			<< rotation(2, 0) << ","
			<< rotation(2, 1) << ","
			<< rotation(2, 2) << ",";
		trajectory_file << center(0) << "," << center(1) << "," << center(2) << "\n";
	}
	trajectory_file.close();
	return true;
}


bool Output_Matchings(std::string filename,const Matches_Provider* matches_provider_) {
  
  std::cout<<"Output_Matchings\n";
  //save the matching after add loop edges
  std::ofstream matching_withloopfile(filename);
  matching_withloopfile << "image_i,image_j,weight\n";
  for (const auto& match_item: matches_provider_->pairWise_matches_)
  {
		matching_withloopfile << match_item.first.first << "," << match_item.first.second << "," << match_item.second.size() << "\n";
  }
  matching_withloopfile.close();
  ////BC end //// 
  return true;
}
  
  /////////////BC start////////////////
bool Output_TriangulatedCorrespondings(std::string filename,const SfM_Data& sfm_data_) {
  
  std::cout<<"Output_TriangulatedCorrespondings\n";
  //save the triangulated correspondings
  std::map<Pair, unsigned int> pair_weights;
  for (const auto& landmark_item : sfm_data_.GetLandmarks())
  {
	  std::set<IndexT> obversed_image;
	  std::transform(landmark_item.second.obs.cbegin(), landmark_item.second.obs.cend()
		  , std::inserter(obversed_image, obversed_image.begin()), stl::RetrieveKey());
	  for (const auto& image_i : obversed_image)
	  {
		  for (const auto& image_j : obversed_image)
		  {
			  if (image_i >= image_j) continue;
			  Pair image_pair(std::min(image_i, image_j), std::max(image_i, image_j));
			  if (pair_weights.count(image_pair)==0)
			  {
				  pair_weights.emplace(image_pair, 0);
			  }
			  pair_weights[image_pair]++;
		  }
	  }
  }
  std::ofstream correspond_file(filename);
  correspond_file << "image_i,image_j,weight\n";
  for (const auto& pair_item : pair_weights)
  {
	  correspond_file << pair_item.first.first << "," << pair_item.first.second << "," << pair_item.second << "\n";
  }
  correspond_file.close();
  return true;
}


bool Output_TriangulatedMatchings(std::string filename,const Matches_Provider* matches_provider_,const SfM_Data& sfm_data_) {
  
  std::cout<<"Output_TriangulatedMatchings\n";
  // triangulated matchings
  std::map<Pair, unsigned int> matching_weights;
  for (const auto& landmark_item : sfm_data_.GetLandmarks())
  {
	  
	  for (const auto& image_i_item : landmark_item.second.obs)
	  {
		  for (const auto& image_j_item : landmark_item.second.obs)
		  {
			  IndexT image_i = image_i_item.first;
			  IndexT image_j = image_j_item.first;
			  if (image_i == image_j) continue;
			  Pair image_pair(image_i,image_j);
			  Pair image_pair_inverse(image_j, image_i);
			  bool find_matching = false;
			  
			  if (matches_provider_->pairWise_matches_.count(image_pair))
			  {
				  for (const matching::IndMatch& match : matches_provider_->pairWise_matches_.at(image_pair))
				  {
					  if ((match.i_ == image_i_item.second.id_feat&&match.j_ == image_j_item.second.id_feat))
					  {
						  find_matching = true;
						  break;
					  }
				  }
			  }
			  else if (matches_provider_->pairWise_matches_.count(image_pair_inverse))
			  {
				  for (const matching::IndMatch& match : matches_provider_->pairWise_matches_.at(image_pair_inverse))
				  {
					  if ((match.j_ == image_i_item.second.id_feat&&match.i_ == image_j_item.second.id_feat))
					  {
						  find_matching = true;
						  break;
					  }
				  }
			  }

			  if (find_matching)
			  {
				  if (matching_weights.count(image_pair) == 0)
				  {
					  matching_weights.emplace(image_pair, 0);
				  }
				  matching_weights[image_pair]++;
			  }
			  
		  }
	  }
  }
  std::ofstream matching_file(filename);
  matching_file << "image_i,image_j,weight\n";
  for (const auto& matching_item : matching_weights)
  {
	  matching_file << matching_item.first.first << "," << matching_item.first.second << "," << matching_item.second << "\n";
  }
  matching_file.close();
  
  

  return true;
}

bool Output_AngleBetweenRotations(std::string filename,const Matches_Provider* matches_provider_,const SfM_Data& sfm_data_) 
{
  std::ofstream file(filename);

  std::vector<double> vec_angularErrors;
  for (const auto& match_item: matches_provider_->pairWise_matches_)
  {
      const std::shared_ptr<View> view_i = sfm_data_.GetViews().at(match_item.first.first);
      const std::shared_ptr<View> view_j = sfm_data_.GetViews().at(match_item.first.second);
	  if (!sfm_data_.GetPoses().count(match_item.first.first) || !sfm_data_.GetPoses().count(match_item.first.second))
	  {
		  if (!sfm_data_.GetPoses().count(match_item.first.first))
		  {
			  std::cout << "View " << match_item.first.first << " is not calibrated\n";
		  }
		  if (!sfm_data_.GetPoses().count(match_item.first.second))
		  {
			  std::cout << "View " << match_item.first.second << " is not calibrated\n";
		  }
		  continue;
	  }
      const Mat3 R_i = sfm_data_.GetPoseOrDie(view_i.get()).rotation();
      const Mat3 R_j = sfm_data_.GetPoseOrDie(view_j.get()).rotation();
      const double angularErrorDegree = R2D(getRotationMagnitude(R_i * R_j.transpose()));
      vec_angularErrors.push_back(angularErrorDegree);
      file<<match_item.first.first<<"-"<<match_item.first.second<<":"<<angularErrorDegree<<"\n";
  }
  std::sort(vec_angularErrors.begin(), vec_angularErrors.end());
  MinMaxMedianMean(file,vec_angularErrors);
  file.close();
}

void MinMaxMedianMean(std::ofstream& infofile,const std::vector<double>& vec_dif)
{
  if(vec_dif.size()==0) 
  {
    infofile<<"Warning:the vector is empty\n";
    return;
  }
  double sum = vec_dif[0];
  for(size_t i = 1 ;i <vec_dif.size(); i++ )
  {
      sum += vec_dif[i];
      if(vec_dif[i]<vec_dif[i-1])
      {
        infofile<<"Error:the difference vector is not sorted\n";
        return;
      }
  }
  infofile<<"Min:"<<vec_dif[0]<<"\n";
  infofile<<"Max:"<<*(vec_dif.rbegin())<<"\n";
  infofile<<"Median:"<<vec_dif[vec_dif.size()/2]<<"\n";
  infofile<<"Mean:"<<sum/((double)vec_dif.size())<<"\n";

}


} // namespace sfm
} // namespace openMVG
