#ifndef TAG_FOURPOINTPOSE_H_
#define TAG_FOURPOINTPOSE_H_

#include <vector>

#include <TooN/TooN.h>
#include <TooN/se3.h>

namespace tag {

/// @defgroup fourpointpose Pose estimation from 4 2D-3D point correspondences
/// This group contains a function and a related RANSAC estimator to compute
/// a camera pose from 4 2D-3D point correspondences.

/// The main function for pose estimation from 4 2D - 3D point correspondences.
/// It implements the algorithm given by
/// Input is a list of 3D positions and a list of 3D vectors describing the pixels on the image plane.
/// Ouput is the SE3 describing the camera pose and a flag signaling if the result is valid.
/// @param[in] points a vector containing 4 3D points
/// @param[in] pixels a vector containing 4 2D pixels as 3D vectors to allow arbitrary image planes
/// @param[out] valid output argument, it is set to true to signal a valid result and false otherwise
/// @return SE3 describing the camera pose
/// @ingroup fourpointpose
TooN::SE3 fourPointPose( const std::vector<TooN::Vector<3> > & points, const std::vector<TooN::Vector<3> > & pixels, bool & valid );

/// A RANSAC estimator using the @ref fourPointPose function. The
/// Correspondence datatype must provide a member position for the 3D point and
/// a member pixel for the 2D pixel.
/// @code
/// struct Correspondence {
///     TooN::Vector<3> position;
///     TooN::Vector<2> pixel;
/// };
/// @endcode
/// @ingroup fourpointpose
/// @ingroup ransac
template <class Correspondence, int ImagePlaneZ = 1>
struct Point4SE3Estimation {
    TooN::SE3 T;
    bool valid;

    Point4SE3Estimation() : valid(false) {  }

    void estimate(const std::vector<Correspondence> & matches) {
        assert(matches.size() >= 4);
        valid = true;

        std::vector<TooN::Vector<3> > points(4);
        std::vector<TooN::Vector<3> > pixels(4);
        for(unsigned int i = 0; i < 4; i ++){
            pixels[i] = unproject(matches[i].pixel);
            pixels[i][2] *= ImagePlaneZ;
            points[i] = matches[i].position;
        }
        T = fourPointPose( points, pixels, valid );
    }

    double getSqError(const Correspondence& m) const {
        if(valid){
            TooN::Vector<3> pos = T * m.position;
            TooN::Vector<2> diff = project(pos) - m.pixel / ImagePlaneZ;
            double disp = diff*diff;
            return disp;
        }
        return 100;
    }
};

} // namespace tag

#endif /*FOURPOINTPOSE*/
