#ifndef TAG_HELPERS_H
#define TAG_HELPERS_H

#include <TooN/se3.h>
#include <TooN/wls.h>
#include <cassert>

namespace tag {

/// @defgroup helpersgroup Helper functions
/// This group contains various helper functions used throughout tag.
/// Most functions come in pairs where one version writes the result to a passed in parameter and
/// the second uses the first to return a new result to the caller.

template <class T> inline const typename T::first_type& first_point(const T& t) { return t.first; }
template <class T> inline const typename T::second_type& second_point(const T& t) { return t.second; }
template <class T> inline double noise(const T& t) { return 1.0; }

/// estimates a 2D homography between two sets of correspondences.
/// The observations passed (via iterators) to the estimate method must allow:
/// @code
/// TooN::Vector<2> a = first_point(*it); // default value is "(*it).first"
/// TooN::Vector<2> b = second_point(*it); // default value is "(*it).second"
/// double R = noise(*it); // default value is "1.0"
/// @endcode
/// At least 4 correspondences are needed. The resulting H returned satisfies b = H a, where a, b
/// are in homogenous coordinates.
/// @param begin iterator pointing to the begin of the sequence
/// @param end iterator to the end of the sequence (after the last element)
/// @param[out] H resulting homography
/// @ingroup helpersgroup

template <class It> void getProjectiveHomography(It begin, It end, TooN::Matrix<3>& H){
    assert(std::distance(begin,end) >= 4);

    TooN::WLS<8> wls;
    for (It it=begin; it!=end; it++) {
        const TooN::Vector<2>& a = first_point(*it);
        const TooN::Vector<2>& b = second_point(*it);
        const TooN::Vector<8> J1 = TooN::makeVector(a[0], a[1], 1, 0, 0, 0, -b[0]*a[0], -b[0]*a[1]);
        const TooN::Vector<8> J2 = TooN::makeVector(0, 0, 0, a[0], a[1], 1, -b[1]*a[0], -b[1]*a[1]);
        wls.add_mJ(b[0], J1, noise(*it));
        wls.add_mJ(b[1], J2, noise(*it));
    }
    wls.compute();
    TooN::Vector<8> h = wls.get_mu();
    H[0] = h.template slice<0,3>();
    H[1] = h.template slice<3,3>();
    H[2][0] = h[6];
    H[2][1] = h[7];
    H[2][2] = 1;
}


/// return version of the 2D homography estimation between two sets of correspondences.
/// The observations passed (via iterators) to the estimate method must allow:
/// @code
/// TooN::Vector<2> a = first_point(*it); // default value is "(*it).first"
/// TooN::Vector<2> b = second_point(*it); // default value is "(*it).second"
/// double R = noise(*it); // default value is "1.0"
/// @endcode
/// At least 4 correspondences are needed. The resulting H returned satisfies b = H a, where a, b
/// are in homogenous coordinates.
/// @param begin iterator pointing to the begin of the sequence
/// @param end iterator to the end of the sequence (after the last element)
/// @return resulting homography
/// @ingroup helpersgroup
template <class It> TooN::Matrix<3> getProjectiveHomography(It begin, It end){
    TooN::Matrix<3> H;
    getProjectiveHomography(begin, end, H);
    return H;
}

/// creates a cross product matrix M from a 3 vector v, such that for all vectors w, the following holds: v ^ w = M * w
/// @param vec the 3 vector input
/// @param[out] result the 3x3 matrix to set to the cross product matrix
/// @ingroup essentialgroup
template<class V, class M > inline void getCrossProductMatrix( const V & vec, M & result ){
    assert(vec.size() == 3);
    assert(result.num_cols() == 3 && result.num_rows() == 3);
    result(0,0) = 0; result(0,1) = -vec[2]; result(0,2) = vec[1];
    result(1,0) = vec[2]; result(1,1) = 0; result(1,2) = -vec[0];
    result(2,0) = -vec[1]; result(2,1) = vec[0]; result(2,2) = 0;
}

/// creates an returns a cross product matrix M from a 3 vector v, such that for all vectors w, the following holds: v ^ w = M * w
/// @param vec the 3 vector input
/// @return the 3x3 matrix to set to the cross product matrix
/// @ingroup essentialgroup
template<class V> inline TooN::Matrix<3> getCrossProductMatrix( const V & vec ){
    TooN::Matrix<3> result;
    getCrossProductMatrix(vec, result);
    return result;
}

/// creates the essential matrix corresponding to a given transformation
/// @param transform the transformation as SE3
/// @param[out] E the 3x3 matrix set to the essential matrix
/// @ingroup essentialgroup
template<class M> inline void getEssentialMatrix(const TooN::SE3<> & transform , M & E){
    //assert(E.num_cols() == 3 && E.num_rows() == 3);
    const TooN::Vector<3> & t = transform.get_translation();
    const TooN::Matrix<3> & r = transform.get_rotation().get_matrix();
    E[0] = t[1] * r[2] - t[2] * r[1];
    E[1] = t[2] * r[0] - t[0] * r[2];
    E[2] = t[0] * r[1] - t[1] * r[0];
}

/// creates and returns the essential matrix corresponding to a given transformation
/// @param transform the transformation as SE3
/// @return the 3x3 matrix set to the essential matrix
/// @ingroup essentialgroup
inline TooN::Matrix<3> getEssentialMatrix(const TooN::SE3<> & transform ){
    TooN::Matrix<3> E;
    getEssentialMatrix(transform, E);
    return E;
}

} // namespace tag

#endif // TAG_HELPERS_H
