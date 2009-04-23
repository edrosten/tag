#ifndef TAG_HANDEYE_H_
#define TAG_HANDEYE_H_

#include <vector>
#include <TooN/se3.h>

namespace tag {

/// @defgroup handeye Hand Eye calibration
/// contains various functions to solve simultaneously for two transformations 
/// a ring of transformations. For example, world - world and sensor - sensor
/// registration for two sensors. The implementation follows roughly Tsai & Lenz 89.
/// It is only a linear method and does not include any nonlinear optimisation.

/// computes the transformation in a circle of transformations:
/// AB -> X -> CD -> ?, so that X is the transformation from B to C. ? is not
/// computed. AB, CD need to contain at least 3 entries to provide 2 independent motions
/// @param[in] AB a vector of measurements for AB
/// @param[in] CD a vector of measurements for CD
/// @return the SE3<> X 
/// @ingroup handeye
TooN::SE3<> computeHandEyeSingle( const std::vector<TooN::SE3<> > & AB, const std::vector<TooN::SE3<> > & CD );

/// computes the pair of transformations that complete the circle:
/// AB -> X -> CD -> Y, so that X is the transformation from B to C and Y the 
/// transformation from D to A. AB, CD need to contain at least 3 entries to provide 2 independent motions
/// @param[in] AB a vector of measurements for AB
/// @param[in] CD a vector of measurements for CD
/// @return pair of SE3<>s, X and Y
/// @ingroup handeye
std::pair<TooN::SE3<>, TooN::SE3<> > computeHandEye( const std::vector<TooN::SE3<> > & AB, const std::vector<TooN::SE3<> > & CD );

/// computes the rotation in a circle of transformations:
/// AB -> X -> CD -> ?, so that X is the rotation from B to C. ? is not
/// computed. AB, CD need to contain at least 3 entries to provide 2 independent motions
/// @param[in] AB a vector of measurements for AB
/// @param[in] CD a vector of measurements for CD
/// @return the SE3<> X 
/// @ingroup handeye
TooN::SO3<>  computeHandEyeSingle( const std::vector<TooN::SO3<> > & AB, const std::vector<TooN::SO3<> > & CD );

/// computes the pair of rotations that complete the circle:
/// AB -> X -> CD -> Y, so that X is the transformation from B to C and Y the 
/// transformation from D to A. AB, CD need to contain at least 3 entries to provide 2 independent motions
/// @param[in] AB a vector of measurements for AB
/// @param[in] CD a vector of measurements for CD
/// @return pair of SO3<> s, X and Y
/// @ingroup handeye
std::pair<TooN::SO3<> , TooN::SO3<> > computeHandEye( const std::vector<TooN::SO3<> > & AB, const std::vector<TooN::SO3<> > & CD );

}

#endif //  TAG_HANDEYE_H_
