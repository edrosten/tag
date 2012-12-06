#ifndef TAG_ABSORIENT_H_
#define TAG_ABSORIENT_H_

#include <vector>
#ifdef WIN32
#include <tuple>
#else
#include <tr1/tuple>
#endif

#include <TooN/TooN.h>
#include <TooN/sim2.h>
#include <TooN/sim3.h>
#include <TooN/SVD.h>

namespace tag {

/// @defgroup absorient Absolute Orientation
/// contains various functions to calculate rotations and rigid transformations
/// between sets of D-dim correspondences. There are special case functions
/// for 2D and 3D points sets that directly return SO2, SO3, SE2 and SE3 objects.

namespace Internal {

template <int D>
std::pair<TooN::Matrix<D>, TooN::DefaultPrecision> computeOrientationScale( const std::vector<TooN::Vector<D> > & a, const std::vector<TooN::Vector<D> > & b ){
	TooN::SizeMismatch<D,D>::test(a.front().size(), b.front().size());
	const int DIM = a.front().size();
	const size_t N = std::min(a.size(), b.size());

	// compute cross correlations
	TooN::Matrix<D> s(DIM,DIM); 
	s = TooN::Zeros;
	for( size_t i = 0; i < N; i++){
		s += b[i].as_col() * a[i].as_row();
	}
	s /= N;

	// SVD of cross correlation matrix
	TooN::SVD<D> svd(s);

	// build S for rotation matrix
	TooN::Vector<D> S(DIM);
	S = TooN::Ones;

	const TooN::DefaultPrecision eps = 1e-8;
	const TooN::DefaultPrecision ds = determinant_gaussian_elimination(s);
	if(ds < -eps){
		S[DIM-1] = -1;
	} else if(ds < eps) { // close to 0 let U * VT decide
		const TooN::DefaultPrecision duv = TooN::determinant_gaussian_elimination(svd.get_U()) 
										* TooN::determinant_gaussian_elimination(svd.get_VT());
		if(duv <  0)
			S[DIM-1] = -1;
	}

	// compute trace(DS)
	TooN::DefaultPrecision scale = 0;
	for(int i = 0; i < DIM; ++i)
		scale += svd.get_diagonal()[i] * S[i];

	return std::make_pair(svd.get_U() * S.as_diagonal() * svd.get_VT(), scale);
}

}

/// computes the rotation between two sets of D-dim points maximizing b * Ta
/// This function returns only the rotation matrix computed after
/// Shinji Umeyama, Least-squares estimation of transformation parameters 
/// between two point patterns. IEEE PAMI, 13(4):376-380, 1991.
/// @param[in] a vector of D-dim points
/// @param[in] b vector of D-dim points
/// @return DxD matrix R describing the rotation such that b = R a
/// @ingroup absorient
template <int D>
inline TooN::Matrix<D> computeRotation( const std::vector<TooN::Vector<D> > & a, const std::vector<TooN::Vector<D> > & b ){
	std::pair<TooN::Matrix<D>, TooN::DefaultPrecision> Rs = Internal::computeOrientationScale( a, b );
	return Rs.first;
}

/// computes the rotation between two sets of 3D points maximizing b * Ta
/// This function returns only the rotation matrix as an SO3<> computed after
/// Shinji Umeyama, Least-squares estimation of transformation parameters 
/// between two point patterns. IEEE PAMI, 13(4):376-380, 1991.
/// @param[in] a vector of 3D points
/// @param[in] b vector of 3D points
/// @return TooN::SO3 containing the rotation such that b = T a
/// @ingroup absorient
inline TooN::SO3<> computeOrientation( const std::vector<TooN::Vector<3> > & a, const std::vector<TooN::Vector<3> > & b ){
	std::pair<TooN::Matrix<3>, TooN::DefaultPrecision> result = Internal::computeOrientationScale( a, b );
	return TooN::SO3<>(result.first);
}

/// computes the rotation between two sets of 2D points maximizing b * Ta
/// This function returns only the rotation matrix as an SO3<> computed after
/// Shinji Umeyama, Least-squares estimation of transformation parameters 
/// between two point patterns. IEEE PAMI, 13(4):376-380, 1991.
/// @param[in] a vector of 2D points
/// @param[in] b vector of 2D points
/// @return TooN::SO3 containing the rotation such that b = T a
/// @ingroup absorient
inline TooN::SO2<> computeOrientation( const std::vector<TooN::Vector<2> > & a, const std::vector<TooN::Vector<2> > & b ){
	std::pair<TooN::Matrix<2>, TooN::DefaultPrecision> result = Internal::computeOrientationScale( a, b );
	return TooN::SO2<>(result.first);
}

/// directly computes a rotation between two pairs of rays in space maximizing b * T a
/// its about 8x faster then using the general computeOrientation for 2 correspondences.
/// @param[in] a1 first input vector
/// @param[in] b1 first output vector
/// @param[in] a2 second input vector
/// @param[in] b2 second output vector
/// @return TooN::SO3 containing the rotation such that b = T a
/// @ingroup absorient
TooN::SO3<> computeOrientation( const TooN::Vector<3> & a1, const TooN::Vector<3> & b1, const TooN::Vector<3> & a2, const TooN::Vector<3> & b2 );

/// computes the rigid transformation between two corresponding D dimensional point sets after 
/// Shinji Umeyama, Least-squares estimation of transformation parameters 
/// between two point patterns. IEEE PAMI, 13(4):376-380, 1991.
/// The result is a rotation matrix R and a translation vector t
/// that map points from vector a to points from vector b : b[i] = R * a[i] + t
/// @param[in] a vector of D-dim points
/// @param[in] b vector of D-dim points
/// @return std::pair containing R and t such that b[i] = R * a[i] + t
/// @ingroup absorient
template <int D>
std::pair<TooN::Matrix<D>, TooN::Vector<D> > computeAbsoluteOrientation( const std::vector<TooN::Vector<D> > & a, const std::vector<TooN::Vector<D> > & b){
	TooN::SizeMismatch<D,D>::test(a.front().size(), b.front().size());
	const int DIM = a.front().size();
	const size_t N = std::min(a.size(), b.size());
	
	if(N == 1){    // quick special case
		TooN::Matrix<D> R(DIM, DIM);
		R = TooN::Identity;
		return std::make_pair(R, b.front() - a.front());
	}

	// compute centroids
	TooN::Vector<D> ma = TooN::Zeros(DIM), mb = TooN::Zeros(DIM);
	for(unsigned i = 0; i < N; ++i){
		ma += a[i];
		mb += b[i];
	}
	ma /= N;
	mb /= N;

	// compute shifted locations
	std::vector<TooN::Vector<D> > ap(N), bp(N);
	for( unsigned i = 0; i < N; ++i){
		ap[i] = a[i] - ma;
		bp[i] = b[i] - ma;
	}
	
	// put resulting transformation together
	std::pair<TooN::Matrix<D>, TooN::DefaultPrecision> Rs = Internal::computeOrientationScale( ap, bp );
	return std::make_pair(Rs.first, mb - Rs.first * ma);
}

/// overload of @ref computeAbsoluteOrientation that computes the rigid transformation between two corresponding 3D point sets 
/// as an SE3 that maps points from vector a to points from vector b : b[i] = SE3 * a[i]
/// @param[in] a vector of 3D points
/// @param[in] b vector of 3D points
/// @return TooN::SE3 T containing the transformation such that b = T a
/// @ingroup absorient
inline TooN::SE3<> computeAbsoluteOrientation( const std::vector<TooN::Vector<3> > & a, const std::vector<TooN::Vector<3> > & b){
	std::pair<TooN::Matrix<3>, TooN::Vector<3> > Rt = computeAbsoluteOrientation<3>( a, b );
	return TooN::SE3<>(Rt.first, Rt.second);
}

/// overload of @ref computeAbsoluteOrientation that computes the rigid transformation between two corresponding 2D point sets 
/// as an SE2 that maps points from vector a to points from vector b : b[i] = SE2 * a[i]
/// @param[in] a vector of 2D points
/// @param[in] b vector of 2D points
/// @return TooN::SE2 T containing the transformation such that b = T a
/// @ingroup absorient
inline TooN::SE2<> computeAbsoluteOrientation( const std::vector<TooN::Vector<2> > & a, const std::vector<TooN::Vector<2> > & b){
	std::pair<TooN::Matrix<2>, TooN::Vector<2> > Rt = computeAbsoluteOrientation<2>( a, b );
	return TooN::SE2<>(Rt.first, Rt.second);
}

/// computes the rigid transformation between two corresponding D dimensional point sets after 
/// Shinji Umeyama, Least-squares estimation of transformation parameters 
/// between two point patterns. IEEE PAMI, 13(4):376-380, 1991.
/// The result is a rotation matrix R, a scale factor s and a translation vector t
/// that map points from vector a to points from vector b : b[i] = s * R * a[i] + t
/// @param[in] a vector of D-dim points
/// @param[in] b vector of D-dim points
/// @return std::tuple containing R, t and s such that b[i] = s * R * a[i] + t
/// @ingroup absorient
template <int D>
std::tr1::tuple<TooN::Matrix<D>, TooN::Vector<D>, TooN::DefaultPrecision > computeSimilarity( const std::vector<TooN::Vector<D> > & a, const std::vector<TooN::Vector<D> > & b){
	TooN::SizeMismatch<D,D>::test(a.front().size(), b.front().size());
	const int DIM = a.front().size();
	const size_t N = std::min(a.size(), b.size());
	
	if(N == 1){    // quick special case
		TooN::Matrix<D> R(DIM, DIM);
		R = TooN::Identity;
		return std::tr1::make_tuple(R, b.front() - a.front(), 1);
	}

	// compute centroids
	TooN::Vector<D> ma = TooN::Zeros(DIM), mb = TooN::Zeros(DIM);
	for(unsigned i = 0; i < N; ++i){
		ma += a[i];
		mb += b[i];
	}
	ma /= N;
	mb /= N;

	// compute shifted locations
	std::vector<TooN::Vector<D> > ap(N), bp(N);
	for( unsigned i = 0; i < N; ++i){
		ap[i] = a[i] - ma;
		bp[i] = b[i] - ma;
	}
	
	// put resulting transformation together
	std::pair<TooN::Matrix<D>, TooN::DefaultPrecision> Rs = Internal::computeOrientationScale( ap, bp );
		// compute scale
	TooN::DefaultPrecision sa = 0;
	for( unsigned int i = 0; i < N; ++i){
		sa += TooN::norm_sq(ap[i]);
	}
	sa /= N;
	const TooN::DefaultPrecision scale = Rs.second / sa;
	return std::tr1::make_tuple(Rs.first, mb - Rs.first * (scale * ma), scale);
}

/// alternative to @ref computeSimilarity that computes the rigid transformation between two corresponding 3D point sets 
/// as an SE3 and a scale S that maps points from vector a to points from vector b : b[i] = SE3 * S * a[i]
/// @param[in] a vector of 3D points
/// @param[in] b vector of 3D points
/// @return a pair consisting of a TooN::SE3 T and a double S containing the transformation such that b = T * S * a
/// @ingroup absorient
inline TooN::SIM3<> computeSimilarity( const std::vector<TooN::Vector<3> > & a, const std::vector<TooN::Vector<3> > & b){
	std::tr1::tuple<TooN::Matrix<3>, TooN::Vector<3>, TooN::DefaultPrecision > Rts = computeSimilarity<3>(a,b);
	return TooN::SIM3<>(TooN::SO3<>(std::tr1::get<0>(Rts)), std::tr1::get<1>(Rts), std::tr1::get<2>(Rts));
}

/// alternative to @ref computeSimilarity that computes the rigid transformation between two corresponding 2D point sets 
/// as an SE2 and a scale S that maps points from vector a to points from vector b : b[i] = SE2 * S * a[i]
/// @param[in] a vector of 2D points
/// @param[in] b vector of 2D points
/// @return a pair consisting of a TooN::SE2 T and a double S containing the transformation such that b = T * S * a
/// @ingroup absorient
inline TooN::SIM2<> computeSimilarity( const std::vector<TooN::Vector<2> > & a, const std::vector<TooN::Vector<2> > & b){
	std::tr1::tuple<TooN::Matrix<2>, TooN::Vector<2>, TooN::DefaultPrecision > Rts = computeSimilarity<2>(a,b);
	return TooN::SIM2<>(TooN::SO2<>(std::tr1::get<0>(Rts)), std::tr1::get<1>(Rts), std::tr1::get<2>(Rts));
}

/// computes the mean rotation of a set of rotations. This is the rotation R such that R^{-1} * R_i is minimal for all R_i.
/// @param[in] r a vector of rotations
/// @return TooN::SO3 mean rotation of all input rotations
/// @ingroup absorient
TooN::SO3<> computeMeanOrientation( const std::vector<TooN::SO3<> > & r);

/// computes a rotation matrix corresponding to a unit quaternion. The quaternion
/// is in the format (q0,qx,qy,qz).
/// @param[in] q a 4-vector containing the coefficients of a unit quaternion as (q0,qx,qy,qz)
/// @return a 3x3 rotation matrix corresponding to the quaternion
/// @ingroup absorient
TooN::Matrix<3> quaternionToMatrix( const TooN::Vector<4> & q );

} // namespace tag

#endif // TAG_ABSORIENT_H_
