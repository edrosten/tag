#ifndef TAG_FIVE_POINT
#define TAG_FIVE_POINT

#include <vector>
#include <utility>
#ifndef WIN32
#include <tr1/array>
#else
#include <array>
#endif
#include <TooN/TooN.h>
#include <TooN/se3.h>

namespace tag {

/// @defgroup essentialgroup Essential Matrix tools
/// various functions dealing with essential matrices including 5 point estimation after Nister, 
/// reconstruction of R,t from essential matrix after Horn, and construction of E from R,t.


/// Computes essential matrices representing the epipolar geometry between correspondences.
/// This function implements Nister's 5 point algorithm. For all of the input points pairs,
/// each out matrix satisfies the condition:
/// \f[ \
///  \vec{\text{second}}\  E \ \vec{\text{first}} = 0. \
///  \f]
/// @param points a array of pairs of directions in 3 space containing correspondences
/// @param initial an inital value for the transformation used as starting point of the optimization
/// @return the optimized transformation
/// @ingroup essentialgroup
/// @ingroup essentialgroup
std::vector<TooN::Matrix<3> > five_point(const std::tr1::array<std::pair<TooN::Vector<3>, TooN::Vector<3> >, 5> & points);

/// reconstructs possible R,t from essential matrix E.
/// The implementation follows the algorithm in 
/// Recovering Baseline and Orientation from 'Essential' Matrix
/// BKP Horn, Jan 1990
/// @param E essential matrix
/// @return vector with 4 SE3s representing the possible transformations
/// @ingroup essentialgroup
std::vector<TooN::SE3<> > se3_from_E( const TooN::Matrix<3> & E );

/// optimizes a transformation representing the epipolar geometry between correspondences.
/// This function minimizes the algebraic error of the epipolar geometry through non-linear optimization of the
/// rotation and direction of the translation.
/// For all of the input points pairs,
/// each out matrix satisfies the condition:
/// \f[ 
/// \vec{\text{second}}\  E \ \vec{\text{first}} = 0. 
/// \f]
/// @param points a vector of pairs of directions in 3 space containing correspondences
/// @param initial an inital value for the transformation used as starting point of the optimization
/// @return the optimized transformation
/// @ingroup essentialgroup
TooN::SE3<> optimize_epipolar(const std::vector<std::pair<TooN::Vector<3>, TooN::Vector<3> > > & points, const TooN::SE3<> & initial);

/// Given an essential matrix \e E and two points \e p and \e q, this
/// functions computes the reprojection errors given by the squared distance from 
/// \e p to the line defined by \f$ E\vec{q} \f$ and the squared distance from 
/// \e q to the line defined by \f$ E^T\vec{p} \f$. If \e E is not an essential matrix
/// then the errors will not be sensible. This function avoids a sqrt compared to
/// tag::essential_reprojection_errors.
///
///@param E \e E: essential matrix
///@param q \e q: right hand (first) point
///@param p \e p: left hand (second) point
///@returns the reprojection errors
/// @ingroup essentialgroup
std::pair<double, double> essential_reprojection_errors_squared(const TooN::Matrix<3>& E, const TooN::Vector<3>&q, const TooN::Vector<3>& p);

/// Given an essential matrix \e E and two points \e p and \e q, this
/// functions computes the reprojection errors given by the signed distance from 
/// \e p to the line defined by \f$ E\vec{q} \f$ and the signed distance from 
/// \e q to the line defined by \f$ E^T\vec{p} \f$. If \e E is not an essential matrix
/// then the errors will not be sensible.
///
///@param E \e E: essential matrix
///@param q \e q: right hand (first) point
///@param p \e p: left hand (second) point
///@returns the reprojection errors
/// @ingroup essentialgroup
std::pair<double, double> essential_reprojection_errors(const TooN::Matrix<3>& E, const TooN::Vector<3>&q, const TooN::Vector<3>& p);
}

#endif // TAG_FIVE_POINT
