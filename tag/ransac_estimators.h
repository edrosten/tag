#ifndef TAG_RANSAC_ESTIMATORS_H
#define TAG_RANSAC_ESTIMATORS_H

#include <vector>
#include <algorithm>
#include <cassert>

#include <TooN/TooN.h>
#include <TooN/helpers.h>
#include <TooN/Cholesky.h>
#include <TooN/wls.h>
#include <TooN/LU.h>
#include <TooN/SVD.h>
#include <TooN/SymEigen.h>
#include <TooN/se3.h>
#include <TooN/wls.h>

#include <tag/helpers.h>
#include <tag/absorient.h>

namespace tag {

#if 0
namespace essential_matrix {

    template <class M> inline int getValidPair(const TooN::Matrix<3>& R1, const TooN::Matrix<3>& R2, const TooN::Vector<2>& e, double z1, const M& m)
    {
	TooN::Vector<2> dm = second_point(m)-e;
	TooN::Vector<3> ha = TooN::unproject(first_point(m));
	TooN::Vector<3> inf1 = R1*ha;
	TooN::Vector<3> inf2 = R2*ha;
	double zp1 = inf1[2];
	double zp2 = inf2[2];
	TooN::Vector<2> pinf1 = project(inf1);
	TooN::Vector<2> pinf2 = project(inf2);
	if (zp1*dm*(pinf1 - e) >= 0) {
	    // R1
	    if (zp1 < 0)
		return z1 <= 0 ? 0 : 1;
	    // check for sign match
	    return  ((pinf1-second_point(m))*dm*z1 >= 0) ? 0 : 1;
	} else {
	    //R2
	    if (zp2 < 0)
		return z1 <= 0 ? 2 : 3;
	    return ((pinf2-second_point(m))*dm*z1 >= 0) ? 2 : 3;
	}
    }

    template <int N, class Accessor> double determinant(const TooN::FixedMatrix<N,N,Accessor>& M)
    {
	TooN::LU<N> lu(M);
	double det = 1;
	for (int i=0; i<N; i++)
	    det *= lu.get_lu()[i][i];
	return det;
    }

/// RANSAC estimator to compute the essential matrix from a set of 2D-2D correspondences
/// The observations passed (via iterators) to the estimate method must allow:
/// @code
/// TooN::Vector<2> a = first_point(*it); // default value is "(*it).first"
/// TooN::Vector<2> b = second_point(*it); // default value is "(*it).second"
/// double R = noise(*it); // default value is "1.0"
/// @endcode
/// The resulting transformation will map from a -> b.
/// @ingroup ransac

    struct EssentialMatrix {
        /// minimal number of correspondences
        static const int hypothesis_size = 8;

	TooN::Matrix<3> E;
	template <class It> bool estimate(It begin, It end) {
	    TooN::Matrix<9> M = TooN::zeros<9,9>();
	    for (It it=begin; it!= end; ++it) {
		const TooN::Vector<2>& a = first_point(*it);
		const TooN::Vector<2>& b = second_point(*it);
		const double factor = 1.0/noise(*it);
		const double m[9] = {b[0]*a[0], b[0]*a[1], b[1]*a[0], b[1]*a[1], b[0], b[1], a[0], a[1], 1};
		for (int j=0; j<9; j++) {
		    for (int k=j; k<9; k++) {
			M[j][k] += m[j]*m[k] * factor;
		    }
		}
	    }
	    for (int j=0; j<9;j++)
		for (int k=j; k<9; k++)
		    M[k][j] = M[j][k];
	    TooN::Matrix<4> M11 = M.template slice<0,0,4,4>();
	    TooN::Matrix<4,5> M12 = M.template slice<0,4,4,5>();
	    TooN::Matrix<5> M22 = M.template slice<4,4,5,5>();
	    TooN::Cholesky<5> chol(M22);
	    if (chol.get_rank() < 5) {
		if (chol.get_rank() < 3)
		    return false;
		TooN::Matrix<3> R;
		// Translation is zero (choose t = [0,0,1])
		tag::getProjectiveHomography(begin, end, R);
		TooN::SO3::coerce(R);
		E[0] = -R[1];
		E[1] = R[0];
		E[2][0] = E[2][1] = E[2][2] = 0;
		return true;
	    }
	    TooN::Matrix<5,4> K = chol.inverse_times(M12.T());
	    TooN::Matrix<4> Q = M11 - M12*K;
	    TooN::SymEigen<4> eigen(Q);
	    TooN::Vector<4> e1 = eigen.get_evectors()[0];
	    TooN::Vector<5> e2 = -(K * e1);
	    E[0][0] = e1[0];
	    E[0][1] = e1[1];
	    E[1][0] = e1[2];
	    E[1][1] = e1[3];
	    E[0][2] = e2[0];
	    E[1][2] = e2[1];
	    E[2] = e2.template slice<2,3>();

	    TooN::SVD<3> svdE(E);
	    const TooN::Vector<3> temp = (TooN::make_Vector, 1, 1, 0);
	    E = svdE.get_U()*TooN::diagmult(temp,svdE.get_VT());
	    return true;
	}

	template <class Match> inline bool isInlier(const Match& m, double r) const {
	    TooN::Vector<3> line = E.template slice<0,0,3,2>()*first_point(m) + E.T()[2];
	    double dot = line.template slice<0,2>() * second_point(m) + line[2];
	    return (dot*dot <= (line[0]*line[0] + line[1]*line[1]) * r*r * noise(m));
	}

    template <class Match> inline double score(const Match& m) const {
	    TooN::Vector<3> line = E.template slice<0,0,3,2>()*first_point(m) + E.T()[2];
	    double dot = line.template slice<0,2>() * second_point(m) + line[2];
	    return dot*dot / (line[0]*line[0] + line[1]*line[1]);
    }


	/// Decompose the essential matrix into four possible SE3s
	/// @param[in] begin beginning iterator for observations
	/// @param[in] end ending iterator for observations
	/// @param[out] group where to store the membership info (0,1,2, or 3) for each observation
	/// @return four pairs of the form {number of votes, SE3}
	template <class It> std::vector<std::pair<size_t, TooN::SE3> > decompose(It begin, It end, std::vector<int>& group)
	{
	    static const double vals[9]={0,-1,0,1,0,0,0,0,1};
	    static const TooN::Matrix<3> Rz(vals);

	    const size_t N = std::distance(begin,end);
	    assert(group.size() >= N);

	    TooN::SVD<3> svdE(E);
	    TooN::Matrix<3> R1 = svdE.get_U()*Rz.T()*svdE.get_VT();
	    TooN::Matrix<3> R2 = svdE.get_U()*Rz*svdE.get_VT();
	    if (essential_matrix::determinant(R1) < 0) {
		R1 = -1*R1;
		R2 = -1*R2;
	    }
	    TooN::Vector<3> t1 = svdE.get_U().T()[2];
	    TooN::Vector<3> t2 = -t1;
	    TooN::Vector<2> epipole = project(t1); // which is the same as project(t2)
	    std::vector<std::pair<size_t, TooN::SE3> > result(4,std::make_pair(0,TooN::SE3()));
	    int i=0;
	    for (It it = begin; it!=end; ++it, ++i) {
		int index = getValidPair(R1, R2, epipole, t1[2], *it);
		result[index].first++;
		group[i] = index;
	    }
	    result[0].second.get_rotation() = result[1].second.get_rotation() = R1;
	    result[2].second.get_rotation() = result[3].second.get_rotation() = R2;
	    result[0].second.get_translation() = result[2].second.get_translation() = t1;
	    result[1].second.get_translation() = result[3].second.get_translation() = t2;
	    return result;
	}

    };
} // close namespace essential_matrix

#endif

// this is deprecated, use 5 point instead
// using essential_matrix::EssentialMatrix;

/// RANSAC estimator to compute an homography from a set of 2D-2D correspondences
/// The observations passed (via iterators) to the estimate method must allow:
/// @code
/// TooN::Vector<2> a = first_point(*it); // default value is "(*it).first"
/// TooN::Vector<2> b = second_point(*it); // default value is "(*it).second"
/// double R = noise(*it); // default value is "1.0"
/// @endcode
/// The resulting transformation will map from a -> b.
/// @ingroup ransac
struct Homography {
    /// homography
    TooN::Matrix<3> H;
    /// minimal number of correspondences
    static const int hypothesis_size = 4;

    Homography() { H = TooN::Identity; }

    template <class It> void estimate(It begin, It end) {
        tag::getProjectiveHomography(begin, end, H);
    }

    template <class M> inline double score(const M& m) const {
        TooN::Vector<3> a = TooN::unproject(first_point(m));
        const TooN::Vector<2> & b = second_point(m);
        const TooN::Vector<2> disp = TooN::project(H * a)  - b;
        return (disp*disp);
    }

    template <class M> inline bool isInlier(const M& m, double r) const {
        return this->score(m) <= r*r * noise(m);
    }
};

/// RANSAC estimator to compute an affine homography from a set of 2D-2D correspondences
/// The observations passed (via iterators) to the estimate method must allow:
/// @code
/// TooN::Vector<2> a = first_point(*it); // default value is "(*it).first"
/// TooN::Vector<2> b = second_point(*it); // default value is "(*it).second"
/// double R = noise(*it); // default value is "1.0"
/// @endcode
/// The resulting transformation will map from a -> b.
/// @ingroup ransac
struct AffineHomography {
    /// the linear part of the affine transformation
    TooN::Matrix<2> A;
    /// the translation part of the affine transformation
    TooN::Vector<2> t;
    /// minimal number of correspondences
    static const int hypothesis_size = 3;

    AffineHomography() : A(TooN::Zeros), t(TooN::Zeros) {}

    template <class It> void estimate(It begin, It end) {
	TooN::WLS<3> wls_x, wls_y;
	wls_x.clear();
	wls_y.clear();
        for (It it = begin; it!= end; ++it) {
	    const TooN::Vector<2>& a = first_point(*it);
	    const TooN::Vector<2>& b = second_point(*it);
	    const double weight = 1.0 / noise(*it);
	    wls_x.add_mJ(b[0], TooN::unproject(a), weight);
	    wls_y.add_mJ(b[1], TooN::unproject(a), weight);
        }
	wls_x.compute();
	wls_y.compute();
	TooN::Vector<3> Atx = wls_x.get_mu();
	TooN::Vector<3> Aty = wls_y.get_mu();
	A[0] = Atx.template slice<0,2>();
	A[1] = Aty.template slice<0,2>();
	t[0] = Atx[2];
	t[1] = Aty[2];
    }

    template <class M> inline double score(const M& m) const {
	const TooN::Vector<2>& a = first_point(m);
	const TooN::Vector<2>& b = second_point(m);
	const TooN::Vector<2> disp = A*a + t - b;
	return (disp*disp);
    }

    template <class M> inline bool isInlier(const M& m, double r) const {
	return this->score(m) <= r*r * noise(m);
    }
};

/// RANSAC estimator to compute a camera rotation between two sets of rays:
/// The rays are specified as 2D vectors on the camera plane. The minimal
/// set are two corresponding pairs of rays.
/// @code
/// TooN::Vector<2> a = first_point(*it); // default value is "(*it).first"
/// TooN::Vector<2> b = second_point(*it); // default value is "(*it).second"
/// double R = noise(*it); // default value is "1.0"
/// @endcode
/// The resulting transformation will map from a -> b.
/// @ingroup ransac
struct CameraRotation {
    /// homography
    TooN::SO3<> rotation;
    /// minimal number of correspondences
    static const int hypothesis_size = 2;

    template <class It> void estimate(It begin, It end) {
        assert(std::distance(begin, end) == hypothesis_size);
        if(TooN::norm_sq(first_point(begin[0])) == first_point(begin[0]) * first_point(begin[1]) ||
           TooN::norm_sq(second_point(begin[0])) == second_point(begin[0]) * second_point(begin[1]))
           return;
        rotation = tag::computeOrientation(TooN::unproject(first_point(begin[0])),
                                      TooN::unproject(second_point(begin[0])),
                                      TooN::unproject(first_point(begin[1])),
                                      TooN::unproject(second_point(begin[1])));
    }

    template <class M> inline double score(const M& m) const {
        const TooN::Vector<2> disp = TooN::project(rotation * TooN::unproject(first_point(m))) - second_point(m);
        return (disp*disp);
    }

    template <class M> inline bool isInlier(const M& m, double r) const {
        return this->score(m) <= r*r * noise(m);
    }
};

/// RANSAC estimator to compute a plane fitting 3 or more points. 
/// The points are Vector<3> or similar. The minimal set are three points, in
/// which case a fast method is used to compute the plane. For more
/// correspondences a method using SVD is employed. Pass in iterators where each
/// element is represents a Vector<3>.
///
/// The estimated plane will be represented by (n,d), where n is a unit vector representing
/// the plane normal and d the distance to the origin. If n == (0 0 0), then no plane
/// was found.
/// @ingroup ransac
struct PlaneFromPoints {
    /// the plane equation coefficients as homogeneous vector with unit normal, or (0,0,0,1)
    TooN::Vector<4> plane;
    /// minimal number of correspondences
    static const int hypothesis_size = 3;

    PlaneFromPoints() : plane(TooN::Zeros) {}

	template <class It> void estimate(It begin, It end){
		assert(std::distance(begin, end) >= 3);
		if( std::distance(begin, end) == 3 ){  // fast special case
			 const TooN::Vector<3> d1 = *(begin+1) - *begin;
			 const TooN::Vector<3> d2 = *(begin+2) - *begin;
			 plane.template slice<0,3>() = d1 ^ d2;
			 TooN::normalize(plane.template slice<0,3>());
			 plane[3] = -(*begin) * plane.template slice<0,3>();
		} else {
			 TooN::Matrix<> design(std::distance(begin, end), 4);
			 for(It p = begin; p != end; ++p)
				design[p-begin] = TooN::unproject(*p);
			TooN::SVD<> s(design);
			plane = s.get_VT()[3];
			const double d = sqrt(plane.template slice<0,3>() * plane.template slice<0,3>());
			if(d > 1e-10){
				plane /= d;
			} else {
				plane = TooN::makeVector( 0, 0, 0, 1);
			}
		}
	}

	template <class M> inline double score(const M & m) const {
		const double d = plane * TooN::unproject(m);
		return d*d;
	}

	template <class M> inline bool isInlier(const M& m, double r) const {
		return this->score(m) <= r*r * noise(m);
	}
};

} // namespace tag

#endif
