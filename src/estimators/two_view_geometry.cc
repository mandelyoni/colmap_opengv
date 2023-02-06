// Copyright (c) 2023, ETH Zurich and UNC Chapel Hill.
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//
//     * Neither the name of ETH Zurich and UNC Chapel Hill nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
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
// Author: Johannes L. Schoenberger (jsch-at-demuc-dot-de)

#include "estimators/two_view_geometry.h"

#include <unordered_set>

#include "base/camera.h"
#include "base/essential_matrix.h"
#include "base/homography_matrix.h"
#include "base/pose.h"
#include "base/projection.h"
#include "base/triangulation.h"
#include "estimators/essential_matrix.h"
#include "estimators/fundamental_matrix.h"
#include "estimators/homography_matrix.h"
#include "estimators/translation_transform.h"
#include "optim/loransac.h"
#include "optim/ransac.h"
#include "util/random.h"

//#include <opengv/relative_pose/methods.hpp>
#include "/home/yoni1/Documents/colmap/src/opengv/include/opengv/relative_pose/CentralRelativeAdapter.hpp"
#include "/home/yoni1/Documents/colmap/src/opengv/include/opengv/sac/Ransac.hpp"
#include "/home/yoni1/Documents/colmap/src/opengv/include/opengv/triangulation/methods.hpp"
//#include <opengv/sac/Lmeds.hpp>
#include "/home/yoni1/Documents/colmap/src/opengv/include/opengv/sac_problems/relative_pose/CentralRelativePoseSacProblem.hpp"

namespace colmap {
namespace {

FeatureMatches ExtractInlierMatches(const FeatureMatches& matches,
                                    const size_t num_inliers,
                                    const std::vector<char>& inlier_mask) 
{
  FeatureMatches inlier_matches(num_inliers);
  size_t j = 0;
  for (size_t i = 0; i < matches.size(); ++i) {
    if (inlier_mask[i]) {
      inlier_matches[j] = matches[i];
      j += 1;
    }
  }

  return inlier_matches;
}

FeatureMatches ExtractOutlierMatches(const FeatureMatches& matches,
                                     const FeatureMatches& inlier_matches) 
{
  CHECK_GE(matches.size(), inlier_matches.size());

  std::unordered_set<std::pair<point2D_t, point2D_t>> inlier_matches_set;
  inlier_matches_set.reserve(inlier_matches.size());
  for (const auto& match : inlier_matches) {
    inlier_matches_set.emplace(match.point2D_idx1, match.point2D_idx2);
  }

  FeatureMatches outlier_matches;
  outlier_matches.reserve(matches.size() - inlier_matches.size());

  for (const auto& match : matches) {
    if (inlier_matches_set.count(
            std::make_pair(match.point2D_idx1, match.point2D_idx2)) == 0) {
      outlier_matches.push_back(match);
    }
  }

  return outlier_matches;
}

inline bool IsImagePointInBoundingBox(const Eigen::Vector2d& point,
                                      const double minx, const double maxx,
                                      const double miny, const double maxy) {
  return point.x() >= minx && point.x() <= maxx && point.y() >= miny &&
         point.y() <= maxy;
}

}  // namespace

void TwoViewGeometry::Invert() {
  F.transposeInPlace();
  E.transposeInPlace();
  H = H.inverse().eval();

  const Eigen::Vector4d orig_qvec = qvec;
  const Eigen::Vector3d orig_tvec = tvec;
  InvertPose(orig_qvec, orig_tvec, &qvec, &tvec);

  for (auto& match : inlier_matches) {
    std::swap(match.point2D_idx1, match.point2D_idx2);
  }
}


// Yoni - called from  TwoViewGeometryVerifier::RUN()
void TwoViewGeometry::Estimate(const Camera& camera1,
                               const std::vector<Eigen::Vector2d>& points1,
                               const Camera& camera2,
                               const std::vector<Eigen::Vector2d>& points2,
                               const FeatureMatches& matches,
                               const Options& options) 
{
  if (options.force_H_use) 
  {
    EstimateHomography(camera1, points1, camera2, points2, matches, options);
  } 
  else if (camera1.HasPriorFocalLength() && camera2.HasPriorFocalLength()) 
  {
    EstimateCalibrated(camera1, points1, camera2, points2, matches, options);
  } 
  else 
  {
    EstimateUncalibrated(camera1, points1, camera2, points2, matches, options);
  }
}


void TwoViewGeometry::EstimateMultiple(
    const Camera& camera1, const std::vector<Eigen::Vector2d>& points1,
    const Camera& camera2, const std::vector<Eigen::Vector2d>& points2,
    const FeatureMatches& matches, const Options& options) 
{
  FeatureMatches remaining_matches = matches;
  std::vector<TwoViewGeometry> two_view_geometries;
  std::cerr << "whyyyyyyyyyyy"  << std::endl;
  while (true) {
    TwoViewGeometry two_view_geometry;
    two_view_geometry.Estimate(camera1, points1, camera2, points2,
                               remaining_matches, options);
    if (two_view_geometry.config == ConfigurationType::DEGENERATE) {
      break;
    }

    if (options.multiple_ignore_watermark) {
      if (two_view_geometry.config != ConfigurationType::WATERMARK) {
        two_view_geometries.push_back(two_view_geometry);
      }
    } else {
      two_view_geometries.push_back(two_view_geometry);
    }

    remaining_matches = ExtractOutlierMatches(remaining_matches,
                                              two_view_geometry.inlier_matches);
  }

  if (two_view_geometries.empty()) {
    config = ConfigurationType::DEGENERATE;
  } else if (two_view_geometries.size() == 1) {
    *this = two_view_geometries[0];
  } else {
    config = ConfigurationType::MULTIPLE;

    for (const auto& two_view_geometry : two_view_geometries) {
      inlier_matches.insert(inlier_matches.end(),
                            two_view_geometry.inlier_matches.begin(),
                            two_view_geometry.inlier_matches.end());
    }
  }
}

bool TwoViewGeometry::EstimateRelativePose(
    const Camera& camera1, const std::vector<Eigen::Vector2d>& points1,
    const Camera& camera2, const std::vector<Eigen::Vector2d>& points2, const bool use_opengv_flag) 
{
  std::cerr << "EstimateRelativePose Entry" << std::endl;
  // We need a valid epipolar geometry to estimate the relative pose.
  if (config != CALIBRATED && config != UNCALIBRATED && config != PLANAR &&
      config != PANORAMIC && config != PLANAR_OR_PANORAMIC) 
  {
    return false;
  }

  opengv::points_t  points;   
  Eigen::Vector2d p_u;
  double theta, phi;
  opengv::bearingVectors_t  bearingVectors1, bearingVectors2;
  Eigen::Vector3d P1, P2; 

  // Extract normalized inlier points.
  std::vector<Eigen::Vector2d> inlier_points1_normalized;
  inlier_points1_normalized.reserve(inlier_matches.size());
  std::vector<Eigen::Vector2d> inlier_points2_normalized;
  inlier_points2_normalized.reserve(inlier_matches.size());
  unsigned long ii = 0;

  for (const auto& match : inlier_matches) 
  {
    const point2D_t idx1 = match.point2D_idx1;
    const point2D_t idx2 = match.point2D_idx2;

    if (use_opengv_flag) //Yoni 
    {
      p_u = camera1.ImageToWorld(points1[idx1]);   // normalized plane coordinates 
      //std::cerr << "use_opengv_flag in EstimateRelativePose" << std::endl;
      theta = atan2(p_u.norm(), 1);
      if (theta < 1e-10)
      {
          phi = 0.0;
      }
      else
      {
          phi = atan2(p_u(1), p_u(0));
      }
      P1[0] = sin(theta) * cos(phi);
      P1[1] = sin(theta) * sin(phi);
      P1[2] = cos(theta);
      //std::cerr << "P1 = " << P1 << std::endl;
      bearingVectors1.push_back(P1);
      bearingVectors1[ii] = bearingVectors1[ii] / bearingVectors1[ii].norm(); // sphere coordinates
      Eigen::Vector2d P1_xy;
      P1_xy[0] = P1[0];
      P1_xy[1] = P1[1];
      inlier_points1_normalized.push_back(P1_xy); //.push_back(P1_xy); 
      
      p_u = camera2.ImageToWorld(points2[idx2]);      
      theta = atan2(p_u.norm(), 1);
      if (theta < 1e-10)
      {
          phi = 0.0;
      }
      else
      {
          phi = atan2(p_u(1), p_u(0));
      }
      P2[0] = sin(theta) * cos(phi);  //x
      P2[1] = sin(theta) * sin(phi);  //y
      P2[2] = cos(theta);             //z 
      bearingVectors2.push_back(P2);
      bearingVectors2[ii] = bearingVectors2[ii] / bearingVectors2[ii].norm(); // sphere coordinates 
      Eigen::Vector2d P2_xy;
      P2_xy[0] = P2[0];
      P2_xy[1] = P2[1];
      inlier_points2_normalized.push_back(P2_xy);   
      ii++;       
    }
    else
    {
      inlier_points1_normalized.push_back(camera1.ImageToWorld(points1[idx1]));
      inlier_points2_normalized.push_back(camera2.ImageToWorld(points2[idx2]));
    }
  }

  Eigen::Matrix3d R;
  std::vector<Eigen::Vector3d> points3D;

  if (config == CALIBRATED || config == UNCALIBRATED) 
  {
    // Try to recover relative pose for calibrated and uncalibrated
    // configurations. In the uncalibrated case, this most likely leads to a
    // ill-defined reconstruction, but sometimes it succeeds anyways after e.g.
    // subsequent bundle-adjustment etc.
    if (use_opengv_flag)
    {
      std::cerr << "  EstimateRelativePose use_openGV = 1" << std::endl;

      opengv::relative_pose::CentralRelativeAdapter adapter(bearingVectors1, bearingVectors2);
      opengv::sac::Ransac<opengv::sac_problems::relative_pose::CentralRelativePoseSacProblem> ransac;
      std::shared_ptr<opengv::sac_problems::relative_pose::CentralRelativePoseSacProblem> relposeproblem_ptr(new opengv::sac_problems::relative_pose::CentralRelativePoseSacProblem(adapter, opengv::sac_problems::relative_pose::CentralRelativePoseSacProblem::EIGHTPT ) );
      // run ransac
      ransac.sac_model_ = relposeproblem_ptr;
      ransac.threshold_ = 2.0*(1.0 - cos(atan(sqrt(2.0)*0.5/800.0)));
      ransac.max_iterations_ = 100000;//
      //ransac.probability_ = 
      ransac.computeModel();

      
    std::cout << "Ransac needed " << ransac.iterations_ << " iterations and ";
    std::cout << "the number of inliers is: " << ransac.inliers_.size();
    std::cout << std::endl << std::endl;
    // std::cout << "the found inliers are: " << std::endl;
    // for(size_t i = 0; i < ransac.inliers_.size(); i++)
    //   std::cout << ransac.inliers_[i] << " ";
    // std::cout << std::endl << std::endl;


      Eigen::Matrix3x4d proj_mat = ransac.model_coefficients_;
      
      R = proj_mat.leftCols(3);
      tvec = proj_mat.col(3);
      int iterations = 100;

      // opengv::translation_t position = tvec; 
      // opengv::rotation_t rotation = R;
      // opengv::relative_pose::CentralRelativeAdapter adapter(
      //   bearingVectors1,
      //   bearingVectors2,
      //   position,
      //   rotation);

      std::cout << "!!!the opengv R results is: " << std::endl;
      std::cout << R << std::endl;
      std::cout << "!!!the normalized opengv translation is: " << std::endl;
      //std::cout << proj_mat.col(3)/proj_mat.col(3).norm() << std::endl;
      std::cout << tvec/tvec.norm() << std::endl;

      std::cout << "running opengv triangulation algorithm 1" << std::endl;
      Eigen::MatrixXd triangulate_results(3, ii);
      for(unsigned long j = 0; j < ii; j++)
      {
        //gettimeofday( &tic, 0 );
        for(int i = 0; i < iterations; i++)
        {
          triangulate_results.block<3,1>(0,j) = opengv::triangulation::triangulate(adapter,j);
          //std::cout << "triangulate_results.block " << j << " = " << triangulate_results.block<3,1>(0,j) << std::endl;  
        }
        // gettimeofday( &toc, 0 );
        // triangulate_time = TIMETODOUBLE(timeval_minus(toc,tic)) / iterations;
      }
      std::cout << "triangulate_results = "  << std::endl; 
      std::cout <<  triangulate_results << std::endl; 
      for (unsigned long j=0; j < ii;  j++)   
      {
        points3D.push_back(triangulate_results.block<3,1>(0,j));
      }
    }
    else
    {
      std::cerr << "  EstimateRelativePose use_openGV = 0 " << std::endl;
      PoseFromEssentialMatrix(E, inlier_points1_normalized, inlier_points2_normalized, &R, &tvec, &points3D);
      std::cout << "!!!the colmap R result is: " << std::endl;
      std::cout << R << std::endl;
      std::cout << "!!!the normalized translation is: " << std::endl;
      std::cout << tvec/tvec.norm() << std::endl ;
      std::cout << std::endl;
    }
  } 
  else if (config == PLANAR || config == PANORAMIC || config == PLANAR_OR_PANORAMIC) 
  {
    Eigen::Vector3d n;
    std::cerr << "WHYYYYYYYYYYYY????????????" << std::endl;
    PoseFromHomographyMatrix(
        H, camera1.CalibrationMatrix(), camera2.CalibrationMatrix(),
        inlier_points1_normalized, inlier_points2_normalized, &R, &tvec, &n,
        &points3D);
  } 
  else 
  {
    return false;
  }


  qvec = RotationMatrixToQuaternion(R); // Yoni123

  if (points3D.empty()) 
  {
    std::cerr << "points3D.empty()" << std::endl;
    tri_angle = 0;
  } 
  else 
  {
    tri_angle = Median(CalculateTriangulationAngles(
        Eigen::Vector3d::Zero(), -R.transpose() * tvec, points3D));
    std::cerr << "This is the median: " << tri_angle  << std::endl;
  }

  if (config == PLANAR_OR_PANORAMIC) 
  {
    if (tvec.norm() == 0) 
    {
      config = PANORAMIC;
      tri_angle = 0;
    } 
    else 
    {
      config = PLANAR;
    }
  }

  return true;
}


void TwoViewGeometry::EstimateCalibrated(
    const Camera& camera1, const std::vector<Eigen::Vector2d>& points1,
    const Camera& camera2, const std::vector<Eigen::Vector2d>& points2,
    const FeatureMatches& matches, const Options& options) 
{
  options.Check();

  if (matches.size() < options.min_num_inliers) 
  {
    config = ConfigurationType::DEGENERATE;
    return;
  }

  std::vector<char> maskk(matches.size());
  for (size_t i = 0; i < matches.size(); ++i)
  {
    maskk[i] = false;
  }

  // Extract corresponding points.
  std::vector<Eigen::Vector2d> matched_points1(matches.size());
  std::vector<Eigen::Vector2d> matched_points2(matches.size());
  std::vector<Eigen::Vector2d> matched_points1_normalized(matches.size());
  std::vector<Eigen::Vector2d> matched_points2_normalized(matches.size());

  opengv::points_t  points;   
  Eigen::Vector2d p_u;
  double theta, phi;
  opengv::bearingVectors_t  bearingVectors1, bearingVectors2;
  Eigen::Vector3d P1, P2;  

  for (size_t i = 0; i < matches.size(); ++i) 
  {  
    const point2D_t idx1 = matches[i].point2D_idx1;
    const point2D_t idx2 = matches[i].point2D_idx2;
    matched_points1[i] = points1[idx1];
    matched_points2[i] = points2[idx2];

    if (options.use_opengv) //Yoni
    {
      p_u = camera1.ImageToWorld(points1[idx1]);   // normalized plane coordinates 
      //std::cerr << "p_u = " << p_u << std::endl;
      theta = atan2(p_u.norm(), 1);
      if (theta < 1e-10)
      {
          phi = 0.0;
      }
      else
      {
          phi = atan2(p_u(1), p_u(0));
      }
      P1[0] = sin(theta) * cos(phi);
      P1[1] = sin(theta) * sin(phi);
      P1[2] = cos(theta);
      //std::cerr << "P1 = " << P1 << std::endl;
      bearingVectors1.push_back(P1);
      bearingVectors1[i] = bearingVectors1[i] / bearingVectors1[i].norm(); // sphere coordinates
      Eigen::Vector2d P1_xy;
      P1_xy[0] = P1[0];
      P1_xy[1] = P1[1];
      matched_points1_normalized[i] = P1_xy; //.push_back(P1_xy); 
      
      p_u = camera2.ImageToWorld(points2[idx2]);      
      theta = atan2(p_u.norm(), 1);
      if (theta < 1e-10)
      {
          phi = 0.0;
      }
      else
      {
          phi = atan2(p_u(1), p_u(0));
      }
      P2[0] = sin(theta) * cos(phi);  //x
      P2[1] = sin(theta) * sin(phi);  //y
      P2[2] = cos(theta);             //z
      bearingVectors2.push_back(P2);
      bearingVectors2[i] = bearingVectors2[i] / bearingVectors2[i].norm(); // sphere coordinates  
      Eigen::Vector2d P2_xy;
      P2_xy[0] = P2[0];
      P2_xy[1] = P2[1];
      matched_points2_normalized[i] = P2_xy;      
    }
    else
    {
      matched_points1_normalized[i] = camera1.ImageToWorld(points1[idx1]);
      matched_points2_normalized[i] = camera2.ImageToWorld(points2[idx2]);
    }
  }
  // std::cerr << "matched_points1" << std::endl;
  // for (size_t i = 0; i < matches.size(); ++i) 
  // {  
  //   std::cerr  << matched_points1[i] << "  " << std::endl;
  // }
  // std::cerr << std::endl;

  if (options.use_opengv)
  {
     const std::vector<char>* best_inlier_mask = nullptr;
    size_t num_inliers = 0;
    std::cerr << "  colmap_opengv  use_openGV = 1" << std::endl;
    // create the central relative adapter
    opengv::relative_pose::CentralRelativeAdapter adapter(bearingVectors1, bearingVectors2);
    // create a RANSAC object
    opengv::sac::Ransac<opengv::sac_problems::relative_pose::CentralRelativePoseSacProblem> ransac;
    // create a CentralRelativePoseSacProblem
    // (set algorithm to STEWENIUS, NISTER, SEVENPT, or EIGHTPT)
    std::shared_ptr<opengv::sac_problems::relative_pose::CentralRelativePoseSacProblem> relposeproblem_ptr(new opengv::sac_problems::relative_pose::CentralRelativePoseSacProblem(adapter, opengv::sac_problems::relative_pose::CentralRelativePoseSacProblem::STEWENIUS ) );
    // run ransac
    ransac.sac_model_ = relposeproblem_ptr;
    ransac.threshold_ = 2.0*(1.0 - cos(atan(sqrt(2.0)*0.5/800.0))) * options.ransac_options.max_error;
    ransac.max_iterations_ = 1000000;
    ransac.computeModel();

    //std::cout << "the ransac threshold is: " << ransac.threshold_ << std::endl;
    std::cout << "the ransac results is: " << std::endl;
    std::cout << ransac.model_coefficients_ <<  std::endl;
    std::cout << "the normalized translation is: " << std::endl;
    std::cout << ransac.model_coefficients_.col(3)/ransac.model_coefficients_.col(3).norm() <<  std::endl;
    std::cout << "Ransac needed " << ransac.iterations_ << " iterations and ";
    // std::cout << ransac_time << " seconds" << std::endl << std::endl;
    std::cout << "the number of inliers is: " << ransac.inliers_.size();
    std::cout << std::endl << std::endl;
    std::cout << "the found inliers are: " << std::endl;
    for(size_t i = 0; i < ransac.inliers_.size(); i++)
      std::cout << ransac.inliers_[i] << " ";
    std::cout << std::endl << std::endl;
 
    if (ransac.inliers_.size() < options.min_num_inliers) 
    {
      config = ConfigurationType::DEGENERATE;
      std::cerr << "ConfigurationType::DEGENERATE opengv::" << std::endl;
      std::cerr << " DEGENERATE inlier_matches.size = "  << inlier_matches.size()<< std::endl;
      return;
    }
    else
    {
      num_inliers = ransac.inliers_.size();
      for (size_t i = 0; i < num_inliers; ++i)
      {
        maskk[ransac.inliers_[i]] = true;
      }
      best_inlier_mask = &maskk; //ransac.inliers_;
      config = ConfigurationType::CALIBRATED;
      std::cerr << "  ConfigurationType::CALIBRATED opengv . num_inliers = " << num_inliers <<std::endl;
    }

    if (best_inlier_mask != nullptr) 
    {
      inlier_matches = ExtractInlierMatches(matches, num_inliers, *best_inlier_mask);
      std::cerr << "best_inlier_mask found for opengv" << std::endl;
    }

  }
  else
  {
    // Estimate epipolar models.
    std::cerr << "  colmap_opengv  use_openGV = 0" << std::endl;

    auto E_ransac_options = options.ransac_options; 
    E_ransac_options.max_error =
        (camera1.ImageToWorldThreshold(options.ransac_options.max_error) +
        camera2.ImageToWorldThreshold(options.ransac_options.max_error)) /  2;

    LORANSAC<EssentialMatrixFivePointEstimator,
            EssentialMatrixFivePointEstimator>  E_ransac(E_ransac_options);
    const auto E_report = E_ransac.Estimate(matched_points1_normalized, matched_points2_normalized);
    E = E_report.model;
   

    if ((!E_report.success ) ||
        (E_report.support.num_inliers < options.min_num_inliers)) 
    {
      config = ConfigurationType::DEGENERATE;
      std::cerr << "ConfigurationType::DEGENERATE" << std::endl;
      return;
    }

    const std::vector<char>* best_inlier_mask = nullptr;
    size_t num_inliers = 0;

    num_inliers = E_report.support.num_inliers;
    best_inlier_mask = &E_report.inlier_mask;
    std::cerr << "  Essen. num_inliers =  " <<  num_inliers  <<std::endl;
    std::cerr << "  best_inlier_mask =  " <<  best_inlier_mask  <<std::endl;
  
    config = ConfigurationType::CALIBRATED;
    std::cerr << "  ConfigurationType::CALIBRATED " << std::endl;
    std::cerr << "   E_report.model: " << std::endl;
    std::cerr <<  E  << std::endl;      

    if (best_inlier_mask != nullptr) 
    {
      inlier_matches = ExtractInlierMatches(matches, num_inliers, *best_inlier_mask);

      std::cerr << "best_inlier_mask found for colmap:" << std::endl;
      
      for (size_t i = 0; i < matches.size(); ++i) 
      {
        if ((*best_inlier_mask)[i]) 
        {
           std::cerr << i << " ";
        }
      }
      std::cerr <<  std::endl;

      // if (options.ransac_options. detect_watermark &&
      //     DetectWatermark(camera1, matched_points1, camera2, matched_points2,
      //                     num_inliers, *best_inlier_mask, options)) 
      // {
      //   config = ConfigurationType::WATERMARK;
      //   std::cerr << "  WATERMARK1 " << std::endl;
      // }
    }
  }
}



void TwoViewGeometry::EstimateUncalibrated(
    const Camera& camera1, const std::vector<Eigen::Vector2d>& points1,
    const Camera& camera2, const std::vector<Eigen::Vector2d>& points2,
    const FeatureMatches& matches, const Options& options) 
{
  std::cerr << " not supposed to be here EstimateUncalibrateddddddddddddddddddddd " << std::endl;
  options.Check();  
    return;  
}

void TwoViewGeometry::EstimateHomography(
    const Camera& camera1, const std::vector<Eigen::Vector2d>& points1,
    const Camera& camera2, const std::vector<Eigen::Vector2d>& points2,
    const FeatureMatches& matches, const Options& options) 
{
  options.Check();
   std::cerr << " Noooooo  EstimateHomography " << std::endl;

    return;
 
}

bool TwoViewGeometry::DetectWatermark(
    const Camera& camera1, const std::vector<Eigen::Vector2d>& points1,
    const Camera& camera2, const std::vector<Eigen::Vector2d>& points2,
    const size_t num_inliers, const std::vector<char>& inlier_mask,
    const Options& options) {
  options.Check();

  std::cerr << " Noooooo  DetectWatermark " << std::endl;

  return true ;
}

}  // namespace colmap
