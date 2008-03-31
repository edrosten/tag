#ifndef THREEPOINTPOSE_H
#define THREEPOINTPOSE_H

#include <TooN/se3.h>
#include <vector>

/// The function for pose estimation from three 2D - 3D point correspondences.
/// It implements the algorithm given by the Fischer and Bolles RANSAC paper, 1980.
/// This function assumes that the three points are in general position (not collinear).
/// Input is an array of 3D cartesian positions and an array of 2D vectors that are the perspective projections of the points.
/// Ouput is up to four poses satisfying the input with positive depths (points in front of the camera).
/// @param[in] x an array containing at least 3 points
/// @param[in] z an array containing the perspective projections of the points given by x in the current pose
/// @param[out] poses the vector onto which any valid poses are appended
/// @return the number of  poses appended to the vector

int three_point_pose(const TooN::Vector<3> x[], const TooN::Vector<2> z[], std::vector<TooN::SE3>& poses);

#endif
