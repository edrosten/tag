#include <tag/absorient.h>

#include <cassert>

#include <TooN/helpers.h>
#include <TooN/determinant.h>

namespace tag {

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
