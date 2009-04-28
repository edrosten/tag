#include <tag/absorient.h>

#include <cassert>

#include <TooN/SymEigen.h>
#include <TooN/helpers.h>

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

TooN::SO3<>  computeOrientation( const std::vector<TooN::Vector<3> > & a, const std::vector<TooN::Vector<3> > & b ){
    const size_t N = a.size();
    // compute cross correlations
    const int x = 0, y = 1, z = 2;
    TooN::Matrix<3> s = TooN::Zeros;
    for( unsigned int i = 0; i < N; i++){
        s += a[i].as_col() * b[i].as_row();
    }

    // create symmetric M for eigenvalue analysis
    TooN::Matrix<4> M;
    M[0][0] = s[x][x] + s[y][y] + s[z][z];
    M[1][0] = M[0][1] = s[y][z] - s[z][y];
    M[2][0] = M[0][2] = s[z][x] - s[x][z];
    M[3][0] = M[0][3] = s[x][y] - s[y][x];
    M[1][1] = s[x][x] - s[y][y] - s[z][z];
    M[2][1] = M[1][2] = s[x][y] + s[y][x];
    M[3][1] = M[1][3] = s[z][x] + s[x][z];
    M[2][2] = - s[x][x] + s[y][y] - s[z][z];
    M[3][2] = M[2][3] = s[y][z] + s[z][y];
    M[3][3] = - s[x][x] - s[y][y] + s[z][z];

    // eigenvalue decomposition to find eigenvector to largest eigenvalue
    TooN::SymEigen<4> ev(M);
    TooN::Vector<4> evals = ev.get_evalues();
    int index = 0;
    for(unsigned int i = index+1; i<4; i++)
        if( evals[i] > evals[index] )
            index = i;
    TooN::Vector<4> evec = ev.get_evectors()[index];
    TooN::SO3<>  result;
    result = quaternionToMatrix(evec);
    return result;
}

// computes the orientation from (e1,e2,e3) -> (a,(a^b)^a,a^b), which means that b the second vector is in the a, b plane
static inline TooN::SO3<>  canonicalOrientation( const TooN::Vector<3> & a, const TooN::Vector<3> & b ){
    TooN::Matrix<3> result;
    result.T()[0] = unit(a);
    result.T()[2] = unit(a ^ b);
    result.T()[1] = result.T()[2] ^ result.T()[0];
    return TooN::SO3<> (result);
}

TooN::SO3<>  computeOrientation( const TooN::Vector<3> & a1, const TooN::Vector<3> & b1, const TooN::Vector<3> & a2, const TooN::Vector<3> & b2 ){
    TooN::SO3<>  r1 = canonicalOrientation( a1, a2 );
    TooN::SO3<>  r2 = canonicalOrientation( b1, b2 );
    const TooN::SO3<>  rAB = r2 * r1.inverse();
    r1 = canonicalOrientation( a2, a1 );
    r2 = canonicalOrientation( b2, b1 );
    const TooN::SO3<>  rBA = r2 * r1.inverse();
    const TooN::SO3<>  diff = rBA * rAB.inverse();
    return TooN::SO3<> ::exp(diff.ln() * 0.5) * rAB;
}

TooN::SE3<>  computeAbsoluteOrientation( const std::vector<TooN::Vector<3> > & a, const std::vector<TooN::Vector<3> > & b){
    // std::assert(a.size() <= b.size());
    const size_t N = a.size();

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
