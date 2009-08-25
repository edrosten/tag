#ifndef TAG_FOURPOINTPOSE_H_
#define TAG_FOURPOINTPOSE_H_

#include <vector>

#include <TooN/TooN.h>
#include <TooN/se3.h>

namespace tag {

/// @defgroup fourpointpose DEPRECATED Pose estimation from 4 2D-3D point correspondences - use 3 point instead
/// This group contains a function and a related RANSAC estimator to compute
/// a camera pose from 4 2D-3D point correspondences.

/// The main function for pose estimation from 4 2D - 3D point correspondences.
/// It implements the algorithm given by
/// This function assumes that the 4 points are in general position (not in a single plane) and
/// can solve the problem for all cases, e.g. points do not have to lie in front of the camera etc.
/// Input is a list of 3D positions and a list of 3D vectors describing the ray under which the transformed points are
/// seen by the origin. It does not assume any special camera setup other than the origin being the center of the
/// camera.
/// Ouput is the SE3 describing the camera pose and a flag signaling if the result is valid.
/// @param[in] points a vector containing 4 3D points
/// @param[in] pixels a vector containing 4 2D pixels as 3D vectors to allow arbitrary image planes
/// @param[out] valid output argument, it is set to true to signal a valid result and false otherwise
/// @param[in] angularError maximal angular error that will be tolerated before no result can be computed. The default value
///            corresponds to 90 deg VOW over 640 pixels.
/// @return SE3 describing the camera pose
/// @ingroup fourpointpose
TooN::SE3<> fourPointPose( const std::vector<TooN::Vector<3> > & points, const std::vector<TooN::Vector<3> > & pixels, bool & valid, const double angularError = 0.14 );

/// A special case of the general @ref fourPointPose function which assumes that points are
/// in front of a given camera plane but now may also lie in the same plane (but not all on one line).
/// The frontal assumption is made for the data, but not encoded in any constraint. Care should be taken
/// that the given 3D vectors for the directions are actually pointing towards the points, not away from them.
/// In the general case, all points in a plane has two solutions, if we assume that all points are
/// in front of the camera, there is only one solution. It also is slighlty more optimized, because
/// it does not need to check as many cases.
/// @param[in] points a vector containing 4 3D points
/// @param[in] pixels a vector containing 4 2D pixels as 3D vectors to allow arbitrary image planes
/// @param[out] valid output argument, it is set to true to signal a valid result and false otherwise
/// @param[in] angularError maximal angular error that will be tolerated before no result can be computed. The default value
///            corresponds to 90 deg VOW over 640 pixels.
/// @return SE3 describing the camera pose
/// @ingroup fourpointpose
TooN::SE3<> fourPointPoseFromCamera( const std::vector<TooN::Vector<3> > & points, const std::vector<TooN::Vector<3> > & pixels, bool & valid, const double angularError = 0.14 );

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
template <int ImagePlaneZ = 1>
struct Point4SE3Estimation {
    /// SE3 describing the transformation from world to camera coordinate frame
    TooN::SE3<> T;
    /// was the estimation valid
    bool valid;
    /// angular error to accept in the 4 point pose estimation
    double angularError;
    /// minimal number of correspondences
    static const int hypothesis_size = 4;

    inline Point4SE3Estimation(double ang = 0.14) : valid(false), angularError(ang) {  }

    template<class It> inline bool estimate(It begin, It end) {
        assert(std::distance(begin,end) >= 4);
        valid = true;

        std::vector<TooN::Vector<3> > points(4);
        std::vector<TooN::Vector<3> > pixels(4);
        unsigned int i = 0;
        for(It match = begin; match != end; match++, i++){
            pixels[i] = unproject(match->pixel);
            pixels[i][2] *= ImagePlaneZ;
            points[i] = match->position;
        }
        T = fourPointPoseFromCamera( points, pixels, valid, angularError );
        return valid;
    }

    template<class Obs, class Tol> inline bool isInlier( const Obs& obs, const Tol& tolerance ) const {
        return score(obs) < tolerance;
    }

    template<class Obs> inline double score(const Obs & obs) const {
        if(valid){
            TooN::Vector<3> pos = T * obs.position;
            TooN::Vector<2> diff = project(pos) - obs.pixel / ImagePlaneZ;
            double disp = diff*diff;
            return disp;
        }
        return 100;
    }

    template<class Obs> inline double getSqError(const Obs & obs) const {
        return score(obs);
    }
};

} // namespace tag

#endif /*FOURPOINTPOSE*/
