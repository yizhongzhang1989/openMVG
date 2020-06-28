// This file is part of OpenMVG_IMU , a branch of OpenMVG
// Author: Bao Chong
// Date:2020/06

#include "openMVG/sfm/sfm_data.hpp"
#include "openMVG/sfm/sfm_data_io.hpp"
#include "openMVG/sfm/sfm_data_transform.hpp"
#include "openMVG/stl/stlMap.hpp"

#include "software/SfM/tools_precisionEvaluationToGt.hpp"
#include "software/SfM/SfMPlyHelper.hpp"


#include "openMVG/geometry/Similarity3.hpp"
#include "openMVG/geometry/Similarity3_Kernel.hpp"
#include "openMVG/robust_estimation/robust_estimator_LMeds.hpp"

#include "software/SfM/import/io_readGTInterface.hpp"
#include "software/SfM/import/io_readGTStrecha.hpp"

#include "third_party/cmdLine/cmdLine.h"
#include "third_party/htmlDoc/htmlDoc.hpp"
#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"

#include <cstdlib>

using namespace openMVG;
using namespace openMVG::sfm;
using namespace openMVG::cameras;

// compute the similarity (all input points)
//                          
//
// @param[in]  imu_data_pos       positions of points in global frame 
// @param[in]  sfm_data_pos       positions of points in local frame
// @param[in]  sbRefine           enable usage of BA
// @param[out] sim                      the similarity computed

// (Takenfrom OpenMVG with modification)

bool ComputeSRTandRefine(const std::vector<Vec3>& imu_data_pos, const std::vector<Vec3>& sfm_data_pos,
	openMVG::geometry::Similarity3& sim, bool sbRefine);

// compute the similarity by ransac(3-points)
//                          
//
// @param[in]  imu_data_pos_const       positions of points in global frame 
// @param[in]  sfm_data_pos_const       positions of points in local frame 
// @param[in]  sbRefine                 enable usage of BA
// @param[out] sim                      the similarity computed

// (Takenfrom OpenMVG with modification)
bool RANSACSRT(const std::vector<Vec3>& imu_data_pos_const, const std::vector<Vec3>& sfm_data_pos_const,
	openMVG::geometry::Similarity3& sim, bool sbRefine);

// evaluate the similarity computed and output evaluation file
//                          
//
// @param[in]  vec_camPosGT             the positions of groundtruth 
// @param[in]  vec_camPosComputed       the positions after aligned
// @param[in]  vec_camRotGT             the camera rotation of groundtruth
// @param[in]  vec_camRotComputed       the camera rotation after aligned
// @param[in]  sim                      the SRT computed
// @param[in]  sOutPath                 the output path of evalution file
// @param[out] _htmlDocStream           evalution file
// (Takenfrom OpenMVG without modification)
void EvaluteResult(
	const std::vector<Vec3> & vec_camPosGT,
	const std::vector<Vec3> & vec_camPosComputed,
	const std::vector<Mat3> & vec_camRotGT,
	const std::vector<Mat3> & vec_camRotComputed,
	const openMVG::geometry::Similarity3& sim,
	const std::string & sOutPath,
	htmlDocument::htmlDocumentStream * _htmlDocStream);

//Note:
//*		the all input sfm datas should be computed from the same matches.
//*		the all input sfm datas and imu data should has the same view informations(length, name of image,key).
//*		the `reconsDirectory` is organized by following format: "recon_path1,recon_path2,recon_path3"
//*		only read VIEWS,INTRINSIACS,EXTRINSICS of imu data

int main(int argc, char **argv)
{
  using namespace std;

  CmdLine cmd;

  std::string
    sGTDirectory,
    reconsDirectory,
    sOutDir = "";
  size_t nreconstructions;

  int sestimate_option = 1;
  bool sbRefine = false;

  cmd.add( make_option('i', sGTDirectory, "") );
  cmd.add( make_option('n', nreconstructions, "") );
  cmd.add( make_option('a', reconsDirectory, "") );
  cmd.add( make_option('o', sOutDir, "outdir") );
  cmd.add( make_option('e', sestimate_option, "similarity_estimate_option"));
  cmd.add( make_option('b', sbRefine, "RefineModel"));

  try {
    if (argc == 1) throw std::string("Invalid command line parameter.");
    cmd.process(argc, argv);
  } catch (const std::string& s) {
    std::cerr << "Usage: " << argv[0] << '\n'
      << "[-i|] imu trajectory path(global frame) .\n"
      << "[-n|] number of reconstruction.\n"
      << "[-a|] the different directories of reconstructions spilted by comma(local frame).\n"
      << "[-o|--outdir] path (where the registration statistics will be saved).\n"
	  << "[-e|similarity_estimate_option] similarity_estimate_option(1: estimated by all poses,0: estimated by three poses)\n"
	  << "[-a|RefineModel] enable usage of BA.\n"
      << std::endl;

    std::cerr << s << std::endl;
    return EXIT_FAILURE;
  }

  if (sOutDir.empty())  {
    std::cerr << "\nIt is an invalid output directory" << std::endl;
    return EXIT_FAILURE;
  }

  if (!stlplus::folder_exists(sOutDir))
    stlplus::folder_create(sOutDir);

  //---------------------------------------
  // Load Parameters and data
  //---------------------------------------

  // 1. Initialize the GT:
  //
  std::cout << "////1. Load the IMU data\n";
  
  SfM_Data Merged_data;
  if (!Load(Merged_data, sGTDirectory, ESfM_Data(VIEWS | INTRINSICS | EXTRINSICS))) {
	  std::cerr << std::endl
		  << "The input Merged_data file \"" << sGTDirectory << "\" cannot be read." << std::endl;
	  return EXIT_FAILURE;
  }

  // 2. Load the scene to align
  // the `reconsDirectory` is organized by following format:
  //	"recon_path1,recon_path2,recon_path3"
  std::cout << "////2. Load the scene to align\n";
  std::vector<SfM_Data> sfm_datas;

  size_t found = 0;
  do
  {
	  size_t found_end = reconsDirectory.find(",", found + 1);

	  std::string recons_path;
	  if (found_end == std::string::npos)
	  {
		  recons_path = reconsDirectory.substr(found);
		  found = found_end;
	  }
	  else
	  {
		  recons_path = reconsDirectory.substr(found, found_end - found);
		  found = found_end + 1;
	  }
	  std::cout << "recons_path:" << recons_path << "\n";
	  SfM_Data sfm_data;
	  if (!Load(sfm_data, recons_path, ESfM_Data(VIEWS | INTRINSICS | EXTRINSICS | STRUCTURE))) {
		  std::cerr << std::endl
			  << "The input SfM_Data file \"" << recons_path << "\" cannot be read." << std::endl;
		  
		  return EXIT_FAILURE;
	  }

	  // Assert that GT and loaded scene have the same image count
	  if (sfm_data.GetViews().size() != Merged_data.GetViews().size())
	  {
		  std::cerr << std::endl
			  << "Dataset views count does not match." << "\n"
			  << "#GT views: " << Merged_data.GetViews().size() << "\n"
			  << "#Scene views: " << sfm_data.GetViews().size() << std::endl;
		  return EXIT_FAILURE;
	  }

	  sfm_datas.emplace_back(std::move(sfm_data));

	  
  } while (found!= std::string::npos);

  
  //Assert that there is none valid scene to align
  if(sfm_datas.empty())
  {
      std::cerr << std::endl
      << "There is none valid scene to align." << "\n";
      return EXIT_FAILURE;
  }


  // 3. Find corresponding camera pose data:
  //
  std::cout << "////3. Find corresponding camera pose data\n";
  
  std::map<IndexT, int> imu_align_count;
  
  Vec3 x_lo(1, 0, 0);
  Vec3 y_lo(0, 1, 0);
  Vec3 z_lo(0, 0, 1);
  std::map<std::pair<IndexT,IndexT>, std::vector<IndexT>> point2d2point3d;
  IndexT trackid_count = 0;  // counter of merged tracks
  for (size_t i = 0; i < sfm_datas.size(); i++)
  {
	  std::vector<Vec3> sfm_data_pos, imu_data_pos;
	  std::vector<Mat3> sfm_data_rot, imu_data_rot;  //used for evaluation
	  std::map<IndexT,IndexT> map_pose_id;
	  std::cout << "//3." << i << " scene is mergeing\n";
	  std::cout << "/find common pose\n";
	  std::cout << "/Pose size:"<<sfm_datas[i].GetPoses().size() << "\n";
	  
	  //collect the poses to be aligned
	  //*	rotations
	  //*	positions: camera center, x, y, z axis of camera frame.
	  for (Views::const_iterator iter = sfm_datas[i].GetViews().begin();
		  iter != sfm_datas[i].GetViews().end();
		  ++iter)
	  {
		  const View * sfm_view = iter->second.get();
		  const View * imu_view = Merged_data.GetViews().at(iter->first).get();
		  
		  //Get the image registered
		  if (sfm_datas[i].GetPoses().count(sfm_view->id_pose) == 0) continue;
		  if (!sfm_datas[i].IsPoseAndIntrinsicDefined(sfm_view)) continue;
		  //assure the image in sfm and imu is same
		  if (sfm_view->s_Img_path != imu_view->s_Img_path)
		  {
			  std::cerr << "the image in imu data and sfm data is not the same.\n";
			  return EXIT_FAILURE;
		  }

		  if (Merged_data.GetPoses().count(imu_view->id_pose) == 0)
		  {
			  std::cerr << "IMU data is inconsistent with sfm data\n";
			  return EXIT_FAILURE;
		  }
		  
		  // when ecountering multiple sfm poses refer to the same imu pose,
		  // only accept the first sfm pose.
		  if (map_pose_id.count(imu_view->id_pose))
		  {
			  std::cout << "duplicate pose\n";
			  continue;
		  }
		  
		  map_pose_id.emplace(imu_view->id_pose , sfm_view->id_pose);

		  const geometry::Pose3& sfm_pose = sfm_datas[i].GetPoses().at(sfm_view->id_pose);
		  const geometry::Pose3& imu_pose = Merged_data.GetPoses().at(imu_view->id_pose);

		  sfm_data_pos.emplace_back(sfm_pose.center());
		  imu_data_pos.emplace_back(imu_pose.center());

		  sfm_data_rot.emplace_back(sfm_pose.rotation());
		  imu_data_rot.emplace_back(imu_pose.rotation());
		  
		  const geometry::Pose3& sfm_pose_inverse = sfm_pose.inverse();
		  const geometry::Pose3& imu_pose_inverse = imu_pose.inverse();

		  //x axis of camera frame 
		  sfm_data_pos.emplace_back(sfm_pose_inverse(x_lo));
		  imu_data_pos.emplace_back(imu_pose_inverse(x_lo));

		  //y axis of camera frame
		  sfm_data_pos.emplace_back(sfm_pose_inverse(y_lo));
		  imu_data_pos.emplace_back(imu_pose_inverse(y_lo));

		  //z axis of camera frame
		  sfm_data_pos.emplace_back(sfm_pose_inverse(z_lo));
		  imu_data_pos.emplace_back(imu_pose_inverse(z_lo));

	  }
  
	  
	  // Visual output of the camera location
	  plyHelper::exportToPly(sfm_data_pos, string(stlplus::folder_append_separator(sOutDir) + "sfm_data_"+std::to_string(i)+".ply").c_str());
	  plyHelper::exportToPly(imu_data_pos, string(stlplus::folder_append_separator(sOutDir) + "imu_data_"+ std::to_string(i) +".ply").c_str());

	  // Evaluation
	  htmlDocument::htmlDocumentStream _htmlDocStream("openMVG Quality evaluation_" + std::to_string(i) + ".");
	  

	  //compute SRT
	  openMVG::geometry::Similarity3 sim;
	  if (sestimate_option) {
		  std::cout << "/compute SRT\n";
		  
		  if (!ComputeSRTandRefine(imu_data_pos, sfm_data_pos, sim,sbRefine))
		  {
			  return EXIT_FAILURE;
		  }
		  
	  }
	  else
	  {
		  std::cout << "/compute SRT by ransac\n";
		  if (!RANSACSRT(imu_data_pos, sfm_data_pos, sim, sbRefine))
		  {
			  return EXIT_FAILURE;
		  }
	  }

	  // Evaluation
	  std::cout << "/Evaluation\n";
	  EvaluteResult(imu_data_pos, sfm_data_pos,
		  imu_data_rot, sfm_data_rot,sim,
		  sOutDir, &_htmlDocStream);

	  ////apply similarity transform
	  std::cout << "/apply similarity transform\n";
	  //apply similarity transform on camera poses
	  for (const auto& id_pair : map_pose_id)
	  {
		  const IndexT sfm_pose_id = id_pair.second;
		  const IndexT imu_pose_id = id_pair.first;
		  if (sfm_datas[i].poses.count(sfm_pose_id) == 0)
		  {
			  std::cerr << "Error:no corresponding pose in sfm\n";
			  return EXIT_FAILURE;
		  }
		  if (Merged_data.poses.count(imu_pose_id) == 0)
		  {
			  std::cerr << "Error:no corresponding pose in Merge data\n";
			  return EXIT_FAILURE;
		  }


		  const geometry::Pose3& sfm_pose = sfm_datas[i].poses.at(sfm_pose_id);
		  geometry::Pose3& imu_pose = Merged_data.poses.at(imu_pose_id);

		  imu_pose = sim(sfm_pose);

		  if (imu_align_count.count(imu_pose_id) == 0)
		  {
			  imu_align_count.emplace(imu_pose_id, 1);
		  }
		  else
		  {
			  imu_align_count[imu_pose_id] += 1;
			  //std::cout << "imu:" << imu_pose_id << " aligned repeatedly\n";
		  }
	  }
	  //apply similarity transform on points
	  for (const auto& landmark_entry : sfm_datas[i].GetLandmarks()) //landmark_entry:map{trackid,{lo_x,obs}}
	  {
		  Landmark landmark;  //can use reference?
		  //Observations obs;
		  landmark.X = sim(landmark_entry.second.X);
		  //landmark.obs = landmark_entry.second.obs;
		  ////assure after merging, one feature corresponds to one 3d points,i.e. accept the first observation.
		  for (const auto& obs_entry : landmark_entry.second.obs)  //obs_entry:map{imageid,{lo_x,featid}}
		  {
			  if (Merged_data.GetViews().count(obs_entry.first) == 0)
			  {
				  std::cout << "The image " << obs_entry.first << "is missing in IMU data\n";
				  continue;
			  }
			  std::pair<IndexT, IndexT> obs_pair(obs_entry.first, obs_entry.second.id_feat);
			  if (point2d2point3d.count(obs_pair) == 0)  //?have another simple and safe way 
			  {
				  point2d2point3d.emplace(obs_pair, std::vector<IndexT>());
			  }

			  point2d2point3d.at(obs_pair).emplace_back(trackid_count);
			  
			  if (point2d2point3d.at(obs_pair).size() > 1)
			  {
				  //std::cout << "id_feat " << obs_entry.second.id_feat << "in image " << obs_entry.first << " has one more 3d points\n";
			  }
			  else if(point2d2point3d.at(obs_pair).size() == 1)
			  {
				  
				  landmark.obs.insert(obs_entry);
			  }
		  }
		  //accept the points can be observed.
		  if (landmark.obs.size() > 0)
		  {
			  Merged_data.structure.emplace(trackid_count, landmark);  //?std::move is right
			  trackid_count++;
		  }
		  
	  }

	  ofstream htmlFileStream( string(stlplus::folder_append_separator(sOutDir) +
		"ExternalCalib_Report_" + std::to_string(i) + ".html").c_str());
	  htmlFileStream << _htmlDocStream.getDoc();
	  std::cout << "//3." << i << " scene is merged successfully\n";
  }

  //4 .save data
  std::cout << "////4 .save data\n";
  if (!Save(
	  Merged_data,
	  stlplus::create_filespec(sOutDir, "merged_data.json").c_str(),
	  ESfM_Data(VIEWS | INTRINSICS|EXTRINSICS|STRUCTURE)))
  {
	  return EXIT_FAILURE;
  }
  if (!Save(
	  Merged_data,
	  stlplus::create_filespec(sOutDir, "merged_data.ply").c_str(),
	  ESfM_Data(VIEWS | INTRINSICS | EXTRINSICS | STRUCTURE)))
  {
	  return EXIT_FAILURE;
  }



  return EXIT_SUCCESS;
}

bool ComputeSRTandRefine(const std::vector<Vec3>& imu_data_pos,const std::vector<Vec3>& sfm_data_pos,
						 openMVG::geometry::Similarity3& sim,bool sbRefine)
{
	
	if (imu_data_pos.size() != sfm_data_pos.size()) {
		std::cerr << "Cannot perform registration, vector sizes are different" << std::endl;
		return false;
	}
	std::vector<Vec3> vec_camPosComputed_T;
	Mat3 R;
	Vec3 t;
	double S;
	////START(Author: BC)++++++++++++++++++++++++++++++++++++++++++++++
	if (sbRefine)
	{
		std::cout << "ComputeSRTandRefine\n";
		Mat x1(3, imu_data_pos.size());
		Mat x2(3, imu_data_pos.size());
		for (size_t i = 0; i < imu_data_pos.size(); ++i) {
			x1.col(i) = sfm_data_pos[i];
			x2.col(i) = imu_data_pos[i];
		}
		// Compute rigid transformation p'i = S R pi + t
		if (!openMVG::geometry::FindRTS(x1, x2, &S, &t, &R))
		{
			std::cerr << std::endl
				<< "Estimate failed" << "\n";
			return false;
		}
	}
	else
	{
		std::cout << "ComputeSRT\n";
		if (!computeSimilarity(imu_data_pos, sfm_data_pos, vec_camPosComputed_T, &S, &R, &t))
		{
			std::cerr << std::endl
				<< "Estimate failed" << "\n";
			return false;
		}
	}
	
	sim = openMVG::geometry::Similarity3(geometry::Pose3(R, -R.transpose()* t / S), S);
	//END(Author: BC)===================================================
	//sim = std::move(openMVG::geometry::Similarity3(geometry::Pose3(R, -R.transpose()* t / S), S)); //?std::move is right
	//if (imu_data_pos.size() != sfm_data_pos.size()) {
	//	std::cerr << "Cannot perform registration, vector sizes are different" << std::endl;
	//	return false;
	//}

	//// Move input point in appropriate container
	//Mat x1(3, imu_data_pos.size());
	//Mat x2(3, imu_data_pos.size());
	//for (size_t i = 0; i < imu_data_pos.size(); ++i) {
	//	x1.col(i) = sfm_data_pos[i];
	//	x2.col(i) = imu_data_pos[i];
	//}
	//// Compute rigid transformation p'i = S R pi + t

	//double S;
	//Vec3 t;
	//Mat3 R;
	//openMVG::geometry::FindRTS(x1, x2, &S, &t, &R);

	//sim=std::move(openMVG::geometry::Similarity3(geometry::Pose3(R, -R.transpose()* t / S), S));
}

bool RANSACSRT(const std::vector<Vec3>& imu_data_pos_const, const std::vector<Vec3>& sfm_data_pos_const,
	openMVG::geometry::Similarity3& sim, bool sbRefine)
{
	std::cout << "RANSACSRT\n";
	assert(imu_data_pos.size() == sfm_data_pos.size());
	////START(Author: BC)++++++++++++++++++++++++++++++++++++++++++++++
	std::vector<Vec3> sfm_data_pos(sfm_data_pos_const.begin(), sfm_data_pos_const.end());
	std::vector<Vec3> imu_data_pos(imu_data_pos_const.begin(), imu_data_pos_const.end());
	//END(Author: BC)===================================================
	double pose_center_robust_fitting_error = 0.0;
	openMVG::geometry::Similarity3 sim_to_center;
	//openMVG::geometry::Similarity3 sim;
	bool b_usable = false;

	const Mat X_SfM_Mat = Eigen::Map<Mat>(sfm_data_pos[0].data(), 3, sfm_data_pos.size());
	const Mat X_IMU_Mat = Eigen::Map<Mat>(imu_data_pos[0].data(), 3, imu_data_pos.size());
	geometry::kernel::Similarity3_Kernel kernel(X_SfM_Mat, X_IMU_Mat);
	const double lmeds_median = openMVG::robust::LeastMedianOfSquares(kernel, &sim);

	
	if (lmeds_median != std::numeric_limits<double>::max())
	{
		b_usable = true; 

		// Display the median residual error once the registration is done
		for (Vec3 & pos : sfm_data_pos) // Transform SfM poses for residual computation
		{
			pos = sim(pos);
		}
		Vec residual = (Eigen::Map<Mat3X>(sfm_data_pos[0].data(), 3, sfm_data_pos.size()) - Eigen::Map<Mat3X>(imu_data_pos[0].data(), 3, imu_data_pos.size())).colwise().norm();
		std::sort(residual.data(), residual.data() + residual.size());
		pose_center_robust_fitting_error = residual(residual.size() / 2);
		std::cout << "RANSAC error:" << pose_center_robust_fitting_error << "\n";
		////START(Author: BC)++++++++++++++++++++++++++++++++++++++++++++++
		if (sbRefine)
		{
			std::cout << "Refine\n";
			double S;
			Vec3 t;
			Mat3 R;
			S = sim.scale_;
			R = sim.pose_.rotation();
			t = sim.pose_.translation();


			std::cout << "None Linear Refine\n";
			openMVG::geometry::Refine_RTS(X_SfM_Mat, X_IMU_Mat, &S, &t, &R);
			sim = std::move(openMVG::geometry::Similarity3(geometry::Pose3(R, -R.transpose()* t / S), S));
		}
		//END(Author: BC)===================================================

		// Apply the found transformation to the SfM Data Scene
		//openMVG::sfm::ApplySimilarity(sim, sfm_data);

		//// Move entire scene to center for better numerical stability
		//Vec3 pose_centroid = Vec3::Zero();
		//for (const auto & pose_it : sfm_data.poses)
		//{
		//	pose_centroid += (pose_it.second.center() / (double)sfm_data.poses.size());
		//}
		//sim_to_center = openMVG::geometry::Similarity3(openMVG::sfm::Pose3(Mat3::Identity(), pose_centroid), 1.0);
		//openMVG::sfm::ApplySimilarity(sim_to_center, sfm_data, true);

	}

	return b_usable;
}

void EvaluteResult(
	const std::vector<Vec3> & vec_camPosGT,
	const std::vector<Vec3> & vec_camPosComputed,
	const std::vector<Mat3> & vec_camRotGT,
	const std::vector<Mat3> & vec_camRotComputed,
	const openMVG::geometry::Similarity3& sim,
	const std::string & sOutPath,
	htmlDocument::htmlDocumentStream * _htmlDocStream
)
{
	// Compute global 3D similarity between the camera position
	std::vector<Vec3> vec_camPosComputed_T(vec_camPosGT.size());
	Mat3 R;
	Vec3 t;
	double S;
	S = sim.scale_;
	R = sim.pose_.rotation();
	t = sim.pose_.translation();


	// Compute statistics and export them
	// -a. distance between camera center
	// -b. angle between rotation matrix

	// -a. distance between camera center
	std::vector<double> vec_residualErrors;
	{
		for (size_t i = 0; i < vec_camPosGT.size(); ++i) {
			const double dResidual = (vec_camPosGT[i] - vec_camPosComputed_T[i]).norm();
			vec_residualErrors.push_back(dResidual);
		}
	}

	// -b. angle between rotation matrix
	std::vector<double> vec_angularErrors;
	{
		std::vector<Mat3>::const_iterator iter1 = vec_camRotGT.begin();
		for (std::vector<Mat3>::const_iterator iter2 = vec_camRotComputed.begin();
			iter2 != vec_camRotComputed.end(); ++iter2, ++iter1) {
			const Mat3 R1 = *iter1; //GT
			const Mat3 R2T = *iter2 * R.transpose(); // Computed

			const double angularErrorDegree = R2D(getRotationMagnitude(R1 * R2T.transpose()));
			vec_angularErrors.push_back(angularErrorDegree);
		}
	}

	// Display residual errors :
	std::cout << "\nBaseline residuals (in GT unit)\n";
	copy(vec_residualErrors.begin(), vec_residualErrors.end(), std::ostream_iterator<double>(std::cout, " , "));
	std::cout << "\nAngular residuals (Degree) \n";
	copy(vec_angularErrors.begin(), vec_angularErrors.end(), std::ostream_iterator<double>(std::cout, " , "));

	std::cout << std::endl << "\nBaseline error statistics : \n ";
	minMaxMeanMedian<double>(vec_residualErrors.begin(), vec_residualErrors.end());
	double minB, maxB, meanB, medianB;
	minMaxMeanMedian<double>(vec_residualErrors.begin(), vec_residualErrors.end(), minB, maxB, meanB, medianB);

	std::cout << std::endl << "\nAngular error statistics : \n ";
	minMaxMeanMedian<double>(vec_angularErrors.begin(), vec_angularErrors.end());
	double minA, maxA, meanA, medianA;
	minMaxMeanMedian<double>(vec_angularErrors.begin(), vec_angularErrors.end(), minA, maxA, meanA, medianA);

	// Export camera position (viewable)
	exportToPly(vec_camPosGT, vec_camPosComputed_T,
		stlplus::create_filespec(sOutPath, "camera_Registered", "ply"));

	exportToPly(vec_camPosGT, vec_camPosComputed,
		stlplus::create_filespec(sOutPath, "camera_original", "ply"));

	//-- Export residual to the HTML report
	{
		using namespace htmlDocument;
		_htmlDocStream->pushInfo("<hr>");
		_htmlDocStream->pushInfo(htmlMarkup("h1", "Compare GT camera position and looking direction."));
		_htmlDocStream->pushInfo(" Display per camera after a 3D similarity estimation:<br>");
		_htmlDocStream->pushInfo("<ul><li>Baseline_Residual -> localization error of camera center to GT (in GT unit),</li>");
		_htmlDocStream->pushInfo("<li>Angular_residuals -> direction error as an angular degree error.</li></ul>");

		std::ostringstream os;
		os << "Baseline_Residual=[";
		std::copy(vec_residualErrors.begin(), vec_residualErrors.end(), std::ostream_iterator<double>(os, " "));
		os << "];";
		_htmlDocStream->pushInfo("<hr>");
		_htmlDocStream->pushInfo(htmlDocument::htmlMarkup("pre", os.str()));

		os.str("");
		os << "mean = " << meanB;
		_htmlDocStream->pushInfo("<hr>");
		_htmlDocStream->pushInfo(htmlDocument::htmlMarkup("pre", os.str()));

		os.str("");
		os << "median = " << medianB;
		_htmlDocStream->pushInfo(htmlDocument::htmlMarkup("pre", os.str()));
		_htmlDocStream->pushInfo("<hr>");

		os.str("");
		os << "Angular_residuals=[";
		std::copy(vec_angularErrors.begin(), vec_angularErrors.end(), std::ostream_iterator<double>(os, " "));
		os << "];";
		_htmlDocStream->pushInfo("<br>");
		_htmlDocStream->pushInfo(htmlDocument::htmlMarkup("pre", os.str()));

		os.str("");
		os << "mean = " << meanA;
		_htmlDocStream->pushInfo("<hr>");
		_htmlDocStream->pushInfo(htmlDocument::htmlMarkup("pre", os.str()));

		os.str("");
		os << "median = " << medianA;
		_htmlDocStream->pushInfo(htmlDocument::htmlMarkup("pre", os.str()));
		_htmlDocStream->pushInfo("<hr>");
	}
}


