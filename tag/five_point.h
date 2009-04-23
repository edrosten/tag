#ifndef TAG_FIVE_POINT
#define TAG_FIVE_POINT

#include <vector>
#include <utility>
#include <tr1/array>
#include <TooN/TooN.h>

namespace tag {

std::vector<TooN::Matrix<3> > five_point(std::tr1::array<std::pair<TooN::Vector<3>, TooN::Vector<3> >, 5> points);

}

#endif // TAG_FIVE_POINT
