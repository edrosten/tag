#include <tag/absorient.h>

#include <cassert>

#include <TooN/SymEigen.h>

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

TooN::SO3 computeOrientation( const std::vector<TooN::Vector<3> > & a, const std::vector<TooN::Vector<3> > & b ){
    const size_t N = a.size();
    // compute cross correlations
    const int x = 0, y = 1, z = 2;
    TooN::Matrix<3> s;
    TooN::Zero(s);
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
    TooN::SO3 result;
    result = quaternionToMatrix(evec);
    return result;
}

TooN::SE3 computeAbsoluteOrientation( const std::vector<TooN::Vector<3> > & a, const std::vector<TooN::Vector<3> > & b){
    //std::assert(a.size() == b.size());
    const size_t N = a.size();

    TooN::Vector<3> ma, mb;
    TooN::Zero(ma);
    TooN::Zero(mb);

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
    TooN::SE3 result;
    result.get_rotation() = computeOrientation( ap, bp );
    result.get_translation() = mb - result.get_rotation() * ma;
    return result;
}

#if 0
// what is this actually doing? should give the same result as computeOrientation
TooN::SO3 computeOrientationFitting( const std::vector<TooN::Vector<3> > & a, const std::vector<TooN::Vector<3> > & b ){
    const unsigned int N = a.size();
    // compute cross correlations
    TooN::Matrix<3> m;
    TooN::Zero(m);
    for( unsigned int i = 0; i < N; i++){
        m += a[i].as_col() * b[i].as_row();
    }

    // compute square root of s
    TooN::Matrix<3> s = m.T() * m;
    TooN::SymEigen<3> ev(s);
    for( unsigned int i = 0; i < 3; i++)
        ev.get_evalues()[i] = sqrt(ev.get_evalues()[i]);

    TooN::SO3 result(ev.backsub(m.T()));
    return result;
}
#endif

TooN::SO3 computeMeanOrientation( const std::vector<TooN::SO3> & r){
    const size_t N = r.size();
    std::vector<TooN::SO3> rt(N);
    TooN::SO3 base = r.front();
    TooN::SO3 baseInv = base.inverse();
    TooN::Vector<3> center;
    Zero(center);
    for(unsigned int i = 0; i < N; i++){
        rt[i] = r[i] * baseInv;
        center += rt[i].ln();
    }
    center /= N;
    TooN::SO3 mean = TooN::SO3::exp(center);
    do {
        Zero(center);
        for(unsigned int i = 0; i < N; i++){
            TooN::SO3 diff = rt[i] * mean.inverse();
            center += diff.ln();
        }
        center /= N;
        mean = TooN::SO3::exp(center) * mean;
    } while(center * center > 1e-12);

    return mean * base;
}

}
