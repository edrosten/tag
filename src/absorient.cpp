#include <tag/absorient.h>

#include <cassert>

#include <TooN/svd.h>
#include <TooN/helpers.h>
#include <TooN/determinant.h>

namespace tag {

static const TooN::DefaultPrecision eps = 1e-8;

TooN::Matrix<3> quaternionToMatrix( const TooN::Vector<4> & q ){
    TooN::Matrix<3> result;
    const int w = 0, x = 1, y = 2, z = 3;
    result(0,0) = q[w]*q[w] + q[x]*q[x] - q[y]*q[y] - q[z]*q[z];
    result(0,1) = 2*(q[x]*q[y] - q[w]*q[z]);
    result(1,0) = 2*(q[x]*q[y] + q[w]*q[z]);
    result(0,2) = 2*(q[x]*q[z] + q[w]*q[y]);
    result(2,0) = 2*(q[x]*q[z] - q[w]*q[y]);
    result(1,1) = q[w]*q[w] - q[x]*q[x] + q[y]*q[y] - q[z]*q[z];
    result(1,2) = 2*(q[y]*q[z] - q[w]*q[x]);
    result(2,1) = 2*(q[y]*q[z] + q[w]*q[x]);
    result(2,2) = q[w]*q[w] - q[x]*q[x] - q[y]*q[y] + q[z]*q[z];
    return result;
}

static std::pair<TooN::SO3<>, TooN::DefaultPrecision> computeOrientationScale( const std::vector<TooN::Vector<3> > & a, const std::vector<TooN::Vector<3> > & b ){
	const size_t N = a.size();
	// compute cross correlations
	TooN::Matrix<3> s = TooN::Zeros;
	for( size_t i = 0; i < N; i++){
		s += b[i].as_col() * a[i].as_row();
	}
	s /= N;

	// SVD of cross correlation matrix
	TooN::SVD<3> svd(s);

	// build S for rotation matrix
	TooN::Matrix<3> S = TooN::Identity;

	const TooN::DefaultPrecision ds = determinant_gaussian_elimination(s);
	if(ds < -eps){
		S(2,2) = -1;
	} else if(ds < eps) { // close to 0 let U * VT decide
		const TooN::DefaultPrecision duv = determinant_gaussian_elimination(svd.get_U()) 
										* determinant_gaussian_elimination(svd.get_VT());
		if(duv <  0)
			S(2,2) = -1;
	}

	// compute trace(DS)
	TooN::DefaultPrecision scale = 0;
	for(int i = 0; i < 3; ++i)
		scale += svd.get_diagonal()[i] * S(i,i);

	return std::make_pair(TooN::SO3<>(svd.get_U() * S * svd.get_VT()), scale);
}

TooN::SO3<> computeOrientation( const std::vector<TooN::Vector<3> > & a, const std::vector<TooN::Vector<3> > & b ){
	std::pair<TooN::SO3<>, TooN::DefaultPrecision> result = computeOrientationScale( a, b );
	return result.first;
}

// computes the orientation from (e1,e2,e3) -> (a,(a^b)^a,a^b), which means that b the second vector is in the a, b plane
static inline TooN::SO3<>  canonicalOrientation( const TooN::Vector<3> & a, const TooN::Vector<3> & b ){
    TooN::Vector<3> n = a ^ b;
    if(norm_sq(n) < 1e-30)
	return TooN::SO3<>();
    TooN::Matrix<3> result;
    result.T()[0] = unit(a);
    result.T()[2] = unit(n);
    result.T()[1] = result.T()[2] ^ result.T()[0];
    return TooN::SO3<> (result);
}

TooN::SO3<> computeOrientation( const TooN::Vector<3> & a1, const TooN::Vector<3> & b1, const TooN::Vector<3> & a2, const TooN::Vector<3> & b2 ){
    TooN::SO3<>  r1 = canonicalOrientation( a1, a2 );
    TooN::SO3<>  r2 = canonicalOrientation( b1, b2 );
    const TooN::SO3<>  rAB = r2 * r1.inverse();
    r1 = canonicalOrientation( a2, a1 );
    r2 = canonicalOrientation( b2, b1 );
    const TooN::SO3<>  rBA = r2 * r1.inverse();
    const TooN::SO3<>  diff = rBA * rAB.inverse();
    return TooN::SO3<> ::exp(diff.ln() * 0.5) * rAB;
}

TooN::SE3<> computeAbsoluteOrientation( const std::vector<TooN::Vector<3> > & a, const std::vector<TooN::Vector<3> > & b){
	assert(a.size() <= b.size());
    const size_t N = a.size();

	if(N == 1){    // quick special case
		return TooN::SE3<>(TooN::SO3<>(), b[0] - a[0]);
	}

    TooN::Vector<3> ma = TooN::Zeros, mb = TooN::Zeros;

    // compute centroids
    for(unsigned int i = 0; i < N; i++){
        ma += a[i];
        mb += b[i];
    }
    ma /= N;
    mb /= N;

    // compute shifted locations
    std::vector<TooN::Vector<3> > ap(N), bp(N);
    for( unsigned int i = 0; i < N; i++){
        ap[i] = a[i] - ma;
        bp[i] = b[i] - ma;
    }

    // put resulting transformation together
    TooN::SE3<>  result;
    result.get_rotation() = computeOrientation( ap, bp );
    result.get_translation() = mb - result.get_rotation() * ma;
    return result;
}

std::pair<TooN::SE3<>, TooN::DefaultPrecision> computeSimilarity( const std::vector<TooN::Vector<3> > & a, const std::vector<TooN::Vector<3> > & b){
	assert(a.size() <= b.size());
	const size_t N = a.size();
	
	if(N == 1){    // quick special case
		return std::make_pair(TooN::SE3<>(TooN::SO3<>(), b[0] - a[0]), 1);
	}
	
	TooN::Vector<3> ma = TooN::Zeros, mb = TooN::Zeros;
	
	// compute centroids
	for(unsigned int i = 0; i < N; ++i){
		ma += a[i];
		mb += b[i];
	}
	ma /= N;
	mb /= N;
	
	// compute shifted locations
	std::vector<TooN::Vector<3> > ap(N), bp(N);
	for( unsigned int i = 0; i < N; ++i){
		ap[i] = a[i] - ma;
		bp[i] = b[i] - ma;
	}
	
	// put resulting transformation together 
	TooN::SE3<>  result;
	std::pair<TooN::SO3<>, TooN::DefaultPrecision> rs = computeOrientationScale( ap, bp );
	result.get_rotation() = rs.first;
	
	std::cout << rs.second << std::endl;
	
	// compute scale
	TooN::DefaultPrecision sa = 0;
	for( unsigned int i = 0; i < N; ++i){
		sa += norm_sq(ap[i]);
	}
	sa /= N;
	const TooN::DefaultPrecision scale = rs.second / sa;
	
	result.get_translation() = mb - result.get_rotation() * (scale * ma);
	return std::make_pair(result, scale);
}

TooN::SO3<>  computeMeanOrientation( const std::vector<TooN::SO3<> > & r){
    const size_t N = r.size();
    std::vector<TooN::SO3<> > rt(N);
    TooN::SO3<>  base = r.front();
    TooN::SO3<>  baseInv = base.inverse();
    TooN::Vector<3> center = TooN::Zeros;
    for(unsigned int i = 0; i < N; i++){
        rt[i] = r[i] * baseInv;
        center += rt[i].ln();
    }
    center /= N;
    TooN::SO3<> mean(center);
    do {
        center = TooN::Zeros;
        for(unsigned int i = 0; i < N; i++){
            TooN::SO3<>  diff = rt[i] * mean.inverse();
            center += diff.ln();
        }
        center /= N;
        mean = TooN::SO3<>::exp(center) * mean;
    } while(center * center > 1e-12);

    return mean * base;
}

}
