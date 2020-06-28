// This file is part of OpenMVG_IMU , a branch of OpenMVG
// Author: Bao Chong
// Date:2020/06

#include "openMVG/matching/indMatch.hpp"
#include "openMVG/matching/indMatch_utils.hpp"
#include "openMVG/matching/svg_matches.hpp"
#include "openMVG/image/image_io.hpp"
#include "openMVG/sfm/pipelines/sfm_features_provider.hpp"
#include "openMVG/sfm/pipelines/sfm_matches_provider.hpp"
#include "openMVG/sfm/sfm_data.hpp"
#include "openMVG/sfm/sfm_data_io.hpp"

#include "third_party/cmdLine/cmdLine.h"
#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"
#include "third_party/progress/progress_display.hpp"

#include "third_party/vectorGraphics/svgDrawer.hpp"

#include <cstdlib>
#include <string>
#include <vector>
#include <fstream>
#include <map>

using namespace openMVG;
using namespace openMVG::features;
using namespace openMVG::matching;
using namespace openMVG::sfm;

// Convert HUE color to RGB
inline float hue2rgb(float p, float q, float t)
{
  if (t < 0) t += 1;
  if (t > 1) t -= 1;
  if (t < 1.f/6.f) return p + (q - p) * 6.f * t;
  if (t < 1.f/2.f) return q;
  if (t < 2.f/3.f) return p + (q - p) * (2.f/3.f - t) * 6.f;
  return p;
}

//
// Converts an HSL color value to RGB. Conversion formula
// adapted from http://en.wikipedia.org/wiki/HSL_color_space.
// Assumes h, s, and l are contained in the set [0, 1] and
// returns r, g, and b in the set [0, 255].
void hslToRgb
(
  float h,
  float s,
  float l,
  uint8_t & r,
  uint8_t & g,
  uint8_t & b
)
{
  if (s == 0)
  {
    r = g = b = static_cast<uint8_t>(l * 255.f); // achromatic
  }
  else
  {
    const float q = l < 0.5f ? l * (1 + s) : l + s - l * s;
    const float p = 2.f * l - q;
    r = static_cast<uint8_t>(hue2rgb(p, q, h + 1.f/3.f) * 255.f);
    g = static_cast<uint8_t>(hue2rgb(p, q, h) * 255.f);
    b = static_cast<uint8_t>(hue2rgb(p, q, h - 1.f/3.f) * 255.f);
  }
}

std::string Matches2SVGString_in_motion
(
  const std::string & left_image_path,
  const std::pair<size_t,size_t> & left_image_size,
  const features::PointFeatures & left_features,
  const std::string & right_image_path,
  const std::pair<size_t,size_t> & right_image_size,
  const features::PointFeatures & right_features,
  const matching::IndMatches & matches,
  const bool b_vertical_display = true,
  const double feature_circle_radius = 4.0,
  const double stroke_size = 2.0
)
{
  const size_t svg_w =
    b_vertical_display ?
    std::max(left_image_size.first, right_image_size.first) :
    left_image_size.first + right_image_size.first;
  const size_t svg_h =
    b_vertical_display ?
    left_image_size.second + right_image_size.second :
    std::max(left_image_size.second, right_image_size.second);
  const size_t svg_offset_x = 0;
    // b_vertical_display ?
    // 0 :
    // left_image_size.first;
  const size_t svg_offset_y = 0;
    // b_vertical_display ?
    // left_image_size.second :
    // 0;

  svg::svgDrawer svgStream(svg_w, svg_h);

  // Draw image side by side
  svgStream.drawImage(left_image_path, left_image_size.first, left_image_size.second);
  //svgStream.drawImage(right_image_path, right_image_size.first, right_image_size.second,
  //  b_vertical_display ? 0 : left_image_size.first,
  //  b_vertical_display ? left_image_size.second : 0);

  std::vector<std::string> colors;
  colors.reserve(matches.size());
  // Draw corresponding matches
  // Perform two loop (it helps to recognize the scene when there is many matches):
  // 1. First draw lines
  for (const auto match_it : matches) {
    // Get back linked features
    const features::PointFeature & L = left_features[match_it.i_];
    const features::PointFeature & R = right_features[match_it.j_];
    { // Compute a flashy colour for the correspondence
      std::ostringstream osCol;
      uint8_t r, g, b;
      hslToRgb( (rand()%360) / 360., 1.0, .5, r, g, b);
      osCol << "rgb(" << (int)r <<',' << (int)g << ',' << (int)b <<")";
      colors.push_back(osCol.str());
    }
    // Draw the line between the corresponding feature positions
    svgStream.drawLine(
      L.x(), L.y(),
      R.x() + svg_offset_x, R.y() + svg_offset_y,
      svg::svgStyle().stroke(colors.back(), stroke_size));
  }
  // 2. Then display features circles
  for (size_t i = 0; i < matches.size(); ++i) {
    // Get back linked features
    const features::PointFeature & L = left_features[matches[i].i_];
    const features::PointFeature & R = right_features[matches[i].j_];
    // Draw the features (circle)
    svgStream.drawCircle(
      L.x(), L.y(), feature_circle_radius,
      svg::svgStyle().stroke(colors[i], stroke_size));
    /*svgStream.drawCircle(
      R.x() + svg_offset_x, R.y() + svg_offset_y, feature_circle_radius,
      svg::svgStyle().stroke(colors[i], stroke_size));*/
  }
  return svgStream.closeSvgFile().str();
}


bool Matches2SVG_in_motion

(
  const std::string & left_image_path,
  const std::pair<size_t,size_t> & left_image_size,
  const features::PointFeatures & left_features,
  const std::string & right_image_path,
  const std::pair<size_t,size_t> & right_image_size,
  const features::PointFeatures & right_features,
  const matching::IndMatches & matches,
  const std::string & svg_filename,
  const bool b_vertical_display = true,
  const double feature_circle_radius = 4.0,
  const double stroke_size = 2.0
)
{
  const std::string svg_content =
    Matches2SVGString_in_motion
    (
      left_image_path,
      left_image_size,
      left_features,
      right_image_path,
      right_image_size,
      right_features,
      matches,
      b_vertical_display,
      feature_circle_radius,
      stroke_size
    );
  // Save the SVG file
  std::ofstream svgFile( svg_filename.c_str() );
  if (svgFile.is_open())
  {
    svgFile << svg_content;
    svgFile.close();
    return true;
  }
  return false;
}

int main(int argc, char ** argv)
{
  CmdLine cmd;

  std::string sSfM_Data_Filename;
  std::string sMatchesDir;
  std::string sMatchFile;
  std::string sOutDir = "";
  bool sOutAscii = false,sOutViews = true;

  cmd.add( make_option('i', sSfM_Data_Filename, "input_file") );
  cmd.add( make_option('d', sMatchesDir, "matchdir") );
  cmd.add( make_option('m', sMatchFile, "matchfile") );
  cmd.add( make_option('o', sOutDir, "outdir") );
  cmd.add(make_option('t', sOutAscii, "outascii"));
  cmd.add(make_option('v', sOutViews, "outviews"));

  try {
      if (argc == 1) throw std::string("Invalid command line parameter.");
      cmd.process(argc, argv);
  } catch (const std::string& s) {
      std::cerr << "Export pairwise matches.\nUsage: " << argv[0] << "\n"
      << "[-i|--input_file file] path to a SfM_Data scene\n"
      << "[-d|--matchdir path]\n"
      << "[-m|--sMatchFile filename]\n"
      << "[-o|--outdir path]\n"
	<< "[-t|--outascii output text file of matches]\n"
		  << "[-v|--outviews output views of matches]\n"
      << std::endl;

      std::cerr << s << std::endl;
      return EXIT_FAILURE;
  }

  if (sOutDir.empty())  {
    std::cerr << "\nIt is an invalid output directory" << std::endl;
    return EXIT_FAILURE;
  }
  if (sMatchesDir.empty()) {
    std::cerr << "\nmatchdir cannot be an empty option" << std::endl;
    return EXIT_FAILURE;
  }
  if (sMatchFile.empty()) {
    std::cerr << "\nmatchfile cannot be an empty option" << std::endl;
    return EXIT_FAILURE;
  }

  //---------------------------------------
  // Read SfM Scene (image view names)
  //---------------------------------------
  SfM_Data sfm_data;
  if (!Load(sfm_data, sSfM_Data_Filename, ESfM_Data(VIEWS|INTRINSICS))) {
    std::cerr << std::endl
      << "The input SfM_Data file \""<< sSfM_Data_Filename << "\" cannot be read." << std::endl;
    return EXIT_FAILURE;
  }

  //---------------------------------------
  // Load SfM Scene regions
  //---------------------------------------
  // Init the regions_type from the image describer file (used for image regions extraction)
  const std::string sImage_describer = stlplus::create_filespec(sMatchesDir, "image_describer", "json");
  std::unique_ptr<Regions> regions_type = Init_region_type_from_file(sImage_describer);
  if (!regions_type)
  {
    std::cerr << "Invalid: "
      << sImage_describer << " regions type file." << std::endl;
    return EXIT_FAILURE;
  }

  // Read the features
  std::shared_ptr<Features_Provider> feats_provider = std::make_shared<Features_Provider>();
  if (!feats_provider->load(sfm_data, sMatchesDir, regions_type)) {
    std::cerr << std::endl
      << "Invalid features." << std::endl;
    return EXIT_FAILURE;
  }
  std::shared_ptr<Matches_Provider> matches_provider = std::make_shared<Matches_Provider>();
  if (!matches_provider->load(sfm_data, sMatchFile)) {
    std::cerr << "\nInvalid matches file." << std::endl;
    return EXIT_FAILURE;
  }

  // ------------
  // For each pair, export the matches
  // ------------

  stlplus::folder_create(sOutDir);
  std::cout << "\n Export pairwise matches" << std::endl;
  const Pair_Set pairs = matches_provider->getPairs();
  C_Progress_display my_progress_bar( pairs.size() );
  for (Pair_Set::const_iterator iter = pairs.begin();
    iter != pairs.end();
    ++iter, ++my_progress_bar)
  {
    const size_t I = iter->first;
    const size_t J = iter->second;

    const View * view_I = sfm_data.GetViews().at(I).get();
    const std::string sView_I= stlplus::create_filespec(sfm_data.s_root_path,
      view_I->s_Img_path);
    const View * view_J = sfm_data.GetViews().at(J).get();
    const std::string sView_J= stlplus::create_filespec(sfm_data.s_root_path,
      view_J->s_Img_path);

    // Get corresponding matches
    const std::vector<IndMatch> & vec_FilteredMatches =
      matches_provider->pairWise_matches_.at(*iter);

    if (!vec_FilteredMatches.empty()) {

      // Draw corresponding features
		if (sOutViews)
		{
			const bool bVertical = false;
			std::ostringstream os;
			os << stlplus::folder_append_separator(sOutDir)
				<< iter->first << "_" << iter->second
				<< "_" << vec_FilteredMatches.size() << "_.svg";
			Matches2SVG_in_motion
			(
				sView_I,
				{ view_I->ui_width, view_I->ui_height },
				feats_provider->getFeatures(view_I->id_view),
				sView_J,
				{ view_J->ui_width, view_J->ui_height },
				feats_provider->getFeatures(view_J->id_view),
				vec_FilteredMatches,
				os.str(),
				bVertical
			);
		}
	  if (sOutAscii) {

       std::stringstream sformatstream;
       sformatstream<< iter->first << "_" << iter->second<< "_" << vec_FilteredMatches.size() << "_.txt";
       std::ofstream asciifile(stlplus::create_filespec(sOutDir,sformatstream.str()));
	   if (!asciifile.is_open())
	   {
		   std::cout << "invalid path:" << stlplus::create_filespec(sOutDir, sformatstream.str()) << "\n";
	   }
	   else
	   {
		   const features::PointFeatures & left_features = feats_provider->getFeatures(view_I->id_view);
		   const features::PointFeatures & right_features = feats_provider->getFeatures(view_J->id_view);
		   for (const auto& indexmatch : vec_FilteredMatches)
		   {
			   const features::PointFeature & L = left_features[indexmatch.i_];
			   const features::PointFeature & R = right_features[indexmatch.j_];
			   asciifile << indexmatch.i_ << " " << indexmatch.j_ << " " << L.x() << " " << L.y() << " " << R.x() << " " << R.x() << "\n";
		   }
	   }
       asciifile.close();

	  }
    }
  }
  return EXIT_SUCCESS;
}




