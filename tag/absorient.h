#ifndef TAG_ABSORIENT_H_
#define TAG_ABSORIENT_H_

#include <vector>

#include <TooN/TooN.h>
#include <TooN/se3.h>

namespace tag {

/// @defgroup absorient Absolute Orientation
/// contains various functions to calculate rotations and rigid transformations
/// between sets of 3D correspondences.

/// computes the rotation between two sets of points maximizing b * Ta
/// this function is part of the absolute orientation algorithm after Horn
/// and is used in @ref computeAbsoluteOrientation function.
/// @param[in] a vector of 3D points
/// @param[in] b vector of 3D points
/// @return TooN::SO3 containing the rotation such that b = T a
/// @ingroup absorient
TooN::SO3 computeOrientation( const std::vector<TooN::Vector<3> > & a, const std::vector<TooN::Vector<3> > & b );

/// directly computes a rotation between two pairs of rays in space maximizing b * T a
/// its about 8x faster then using the general computeOrientation for 2 correspondences.
/// @param[in] a1 first input vector
/// @param[in] b1 first output vector
/// @param[in] a2 second input vector
/// @param[in] b2 second output vector
/// @return TooN::SO3 containing the rotation such that b = T a
/// @ingroup absorient
TooN::SO3 computeOrientation( const TooN::Vector<3> & a1, const TooN::Vector<3> & b1, const TooN::Vector<3> & a2, const TooN::Vector<3> & b2 );

/// computes the rigid transformation between two corresponding point sets after Horn
/// The result is an SE3 that maps points from vector a to points from vector b : b[i] = SE3 * a[i]
/// @param[in] a vector of 3D points
/// @param[in] b vector of 3D points
/// @return TooN::SE3 containing the transformation such that b = T a
/// @ingroup absorient
TooN::SE3 computeAbsoluteOrientation( const std::vector<TooN::Vector<3> > & a, const std::vector<TooN::Vector<3> > & b);

/// computes the mean rotation of a set of rotations. This is the rotation R such that R^{-1} * R_i is minimal for all R_i.
/// @param[in] r a vector of rotations
/// @return TooN::SO3 mean rotation of all input rotations
/// @ingroup absorient
TooN::SO3 computeMeanOrientation( const std::vector<TooN::SO3> & r);

/// computes a rotation matrix corresponding to a unit quaternion. The quaternion
/// is in the format (q0,qx,qy,qz) to fit the absolute orientation algorithm. This
/// is a helper function for the @ref computeOrientation function.
/// @param[in] q a 4-vector containing the coefficients of a unit quaternion as (q0,qx,qy,qz)
/// @return a 3x3 rotation matrix corresponding to the quaternion
/// @ingroup absorient
TooN::Matrix<3> quaternionToMatrix( const TooN::Vector<4> & q );

} // namespace tag

#endif // TAG_ABSORIENT_H_
