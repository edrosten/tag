#include <tag/handeye.h>
#include <tag/absorient.h>

#include <TooN/Cholesky.h>

#include <cassert>
#include <iostream>

using namespace std;
using namespace TooN;

namespace tag {

static inline Vector<3> getRotationVector( const SO3<>  & r ){
    const Matrix<3> & my_matrix = r.get_matrix();
    Vector<3> result;
    result[0] = (my_matrix[2][1]-my_matrix[1][2]);
    result[1] = (my_matrix[0][2]-my_matrix[2][0]);
    result[2] = (my_matrix[1][0]-my_matrix[0][1]);
    result = unit(result);
    return result;
}

static inline Vector<3> getRotationVector( const SE3<>  & t ){
    return getRotationVector(t.get_rotation());
}

template<class T>
static inline TooN::SO3<>  solveXABX( const std::vector<T> & A, const std::vector<T> & B ){
    vector<Vector<3> > va(A.size()), vb(A.size());
    for(unsigned int i = 0; i < A.size(); ++i){
        va[i] = getRotationVector(A[i]);
        vb[i] = getRotationVector(B[i]);
    }
    return computeOrientation(va, vb);
}

static inline Matrix<3> eyeMinus( const Matrix<3> & m ){
    Matrix<3> result = Identity;
    result -= m;
    return result;
}

SE3<>  computeHandEyeSingle( const vector<SE3<> > & AB, const vector<SE3<> > & CD ){
    vector<SE3<> > A(AB.size()-1),B(AB.size()-1);
    for(unsigned int i = 0; i < AB.size() - 1; ++i){
        A[i] = CD[i] * CD[i+1].inverse();
        B[i] = AB[i].inverse() * AB[i+1];
    }
    SO3<> R = solveXABX(A,B);
    Matrix<3> JTJ = Zeros;
    Vector<3> JTE = Zeros;
    for(unsigned int i = 0; i < A.size(); ++i){
        Matrix<3> m =  eyeMinus( B[i].get_rotation().get_matrix());
        JTJ += m.T() * m;
        JTE += m.T() * (B[i].get_translation() - R * A[i].get_translation());
    }
    Cholesky<3> chol(JTJ);
    return SE3<> (R, chol.backsub(JTE));
}

std::pair<TooN::SE3<> , TooN::SE3<> > computeHandEye( const std::vector<TooN::SE3<> > & AB, const std::vector<TooN::SE3<> > & CD ){
    assert( AB.size() == CD.size() && AB.size() > 2 );
    return make_pair(computeHandEyeSingle(AB, CD), computeHandEyeSingle(CD, AB)); 
}

SO3<>  computeHandEyeSingle( const vector<SO3<> > & AB, const vector<SO3<> > & CD ){
    vector<SO3<> > A,B;
    for(unsigned int i = 0; i < AB.size() - 1; ++i){
		for(unsigned int j = i+1; j < AB.size(); ++j){
			A.push_back(CD[i] * CD[j].inverse());
			B.push_back(AB[i].inverse() * AB[j]);
		}
    }
    return solveXABX(A,B);
}

std::pair<TooN::SO3<> , TooN::SO3<> > computeHandEye( const std::vector<TooN::SO3<> > & AB, const std::vector<TooN::SO3<> > & CD ){
    assert( AB.size() == CD.size() && AB.size() > 2 );
    return make_pair(computeHandEyeSingle(AB, CD), computeHandEyeSingle(CD, AB));
}

}
