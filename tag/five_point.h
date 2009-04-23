#ifndef TAG_FIVE_POINT
#define TAG_FIVE_POINT

#include <vector>
#include <utility>
#include <tr1/array>
#include <TooN/TooN.h>
#include <TooN/se3.h>

namespace tag {

/// @defgroup essentialgroup Essential Matrix tools
/// various functions dealing with essential matrices including 5 point estimation after Nister, 
/// reconstruction of R,t from essential matrix after Horn, and construction of E from R,t.


/// @ingroup essentialgroup
std::vector<TooN::Matrix<3> > five_point(std::tr1::array<std::pair<TooN::Vector<3>, TooN::Vector<3> >, 5> points);

/// reconstructs possible R,t from essential matrix E.
/// The implementation follows the algorithm in 
/// Recovering Baseline and Orientation from 'Esssential' Matrix
/// BKP Horn, Jan 1990
/// @arg E essential matrix
/// @return vector with 4 SE3s representing the possible transformations
/// @ingroup essentialgroup
std::vector<TooN::SE3<> > se3_from_E( const TooN::Matrix<3> & E );

}

#endif // TAG_FIVE_POINT
