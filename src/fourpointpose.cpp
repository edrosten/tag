#include <tag/fourpointpose.h>
#include <tag/absorient.h>

#include <TooN/SVD.h>
#include <TooN/helpers.h>

#include <iostream>

namespace tag {

// simple helper function to get the two possible ray lengths given a known point, and the distance to the other point
// does various checks to account for numerical inaccuracies
static inline bool computeDistances( const double xx, const double x, const double angle, const double distance, const double angularError, TooN::Vector<2> & roots){
    double p = -x*angle;
    double q = xx - distance;
    double det = p*p*0.25 - q;
    if( det < 0 ){
        double alpha = std::acos(0.5 * angle);
        double beta = std::asin(sqrt(distance) / x);
        double diff = fabs(alpha) - fabs(beta);
        if( diff > angularError ){
			roots = TooN::Zeros;
            return false;
        }
        det = 0;
    } else {
        det = sqrt(det);
    }
    roots = TooN::makeVector( -p*0.5 + det, -p*0.5 - det);
    return true;
}

// helper function for @ref fourPointPose computing the coefficients of the A matrix
static inline TooN::Vector<5> getACoeffs(double c1, double c2, double c3, double d1, double d2, double d3){
    TooN::Vector<5> coeffs;
    // x^4 coefficient
    coeffs[4] = (c3*c3 - c1*c2*c3 - 4 + c1*c1 + c2*c2);
    coeffs[4] = coeffs[4]*coeffs[4];
    // x^3 coefficient
    coeffs[3] = -(c3*c3 - c1*c2*c3 - 4 + c1*c1 + c2*c2)*
                    (-2*d3*c1*c1 + 2*d2*c1*c1 + c2*c2*d3*c1*c1 - c2*c3*c1*(d1 + d2 + d3) +
                    8*(d3 - d2 - d1) + 2*d2*c3*c3 + 2*c3*c3*d1 - 2*d3*c2*c2 + 2*d1*c2*c2 );
    // x^2 coefficient
    coeffs[2] =
        -2*d2*d2*c2*c2 + 24*(d3*d3+d1*d1+d2*d2) + 48*(d2*d1 - d3*d1 - d2*d3) - 28*d2*c3*c3*d1  - 10*d1*d1*c3*c3 + c2*c2*c2*c2*d3*d3
        + 12*d3*d2*c2*c2 - 4*d3*d2*c2*c3*c1 - 2*d1*d1*c1*c1 - 4*d3*c1*c1*d2*c2*c2 + 2*d3*c1*c1*c1*d2*c2*c3 - 2*d3*c1*c1*c1*c1*d2 - 4*c3*c3*c3*d1*d2*c2*c1 - d2*d2*c3*c3*c3*c2*c1
        + d2*d2*c3*c3*c2*c2 + 3*d2*d2*c3*c3*c1*c1 - 2*c2*c2*c2*c2*d3*d1 - 10*(d3*d3*c2*c2 + d1*d1*c2*c2 + d3*d3*c1*c1) -2*d3*d3*c3*c3 + d3*d3*c1*c1*c1*c1
        - 10*d2*d2*c3*c3 - 10*c1*c1*d2*d2 + c3*c3*c3*c3*(d2*d2 + d1*d1) +c2*c2*c2*c2*d1*d1 + d2*d2*c1*c1*c1*c1 + 2*d2*d2*c2*c3*c1 + 12*d1*d3*c1*c1 + 12*d3*d1*c3*c3
        + 20*d3*d1*c2*c2 + 20*d2*d3*c1*c1 + 12*d2*c3*c3*d3 - 4*d2*c3*c3*d3*c1*c1 + 4*d2*c3*c3*c3*c3*d1 + 2*d2*c3*c3*d1*c1*c1 + 2*d2*c3*c3*d1*c2*c2 +2*c2*c3*c1*d3*d3
        - 4*c2*c3*c1*d3*d1 + 2*c2*c3*c1*d1*d1 - 12*d1*d2*c2*c2 +20*d1*d2*c2*c3*c1 - 12*d1*d2*c1*c1 -c3*c2*c1*c1*c1*d3*d3 - c2*c2*c2*c3*c1*d3*d3 - c2*c3*c1*c1*c1*d1*d3
        + c2*c2*c3*c3*d3*d3 - 4*c2*c2*c3*c3*d3*d1 + 3*c2*c2*c3*c3*d1*d1 +3*c2*c2*c1*c1*d3*d3 - 4*c2*c2*c1*c1*d3*d1 + c2*c2*c1*c1*d1*d1 +c2*c2*c3*c3*d1*d3*c1*c1 + 2*c2*c2*c2*c3*c1*d3*d1
        - c2*c2*c2*c3*c1*d1*d1 - d2*d2*c1*c1*c1*c2*c3 + d2*d2*c1*c1*c2*c2 + c2*c2*c1*c1*d2*c3*c3*d1 - c2*c2*c2*d3*d2*c3*c1 - c2*c2*c2*d1*d2*c3*c1 - c2*c3*c3*c3*c1*d3*d2 - c2*c3*c3*c3*c1*d3*d1
        + c2*c2*c3*c3*c1*c1*d3*d2 - c2*c3*c3*c3*c1*d1*d1 - c2*c3*c1*c1*c1*d1*d2 + c3*c3*c1*c1*d3*d3 + c3*c3*c1*c1*d1*d1;
    // x^1 coefficient
    coeffs[1] =
        - 2*d3*d3*d3*c1*c1 - 8*d2*d2*d2 + 14*d2*c3*c3*d1*d1 + 24*(-d1*d1*d2 -d3*d3*d1 - d1*d2*d2 - d3*d3*d2 + d3*d2*d2 + d3*d1*d1)
        - 2*c2*c2*d3*d3*d3 + 2*d2*d2*d2*c1*c1 + 2*c3*c3*d1*d1*d1 + 2*c2*c2*d1*d1*d1 + 2*d2*d2*d2*c3*c3 + 2*c3*c3*d1*d3*d3 - 4*c3*c3*d1*d1*d3 + 2*d2*c3*c3*d3*d3
        - 16*d2*c3*c3*d3*d1 -5*d1*d1*d2*c2*c3*c1 + 2*d1*d1*d2*c1*c1 + c2*c3*c3*c3*c1*d3*d2*d1 + c2*c3*c3*c3*c1*d1*d1*d2
        - d2*c3*c3*c1*c1*d3*d3 +2*d2*d2*c3*c3*c1*c1*d3 - d2*c3*c3*c1*c1*d1*d1
        - d2*d2*d2*c3*c3*c1*c1 + d2*d2*d2*c2*c3*c1 + 2*d2*d2*d1*c2*c2 -5*d2*d2*d1*c2*c3*c1 + 4*d2*d2*d1*c1*c1
        - c2*c2*c3*c3*d1*d3*d3 + 2*c2*c2*c3*c3*d1*d1*d3 - c2*c2*c3*c3*d1*d1*d1 + 6*c2*c2*d3*d3*d1
        - 6*c2*c2*d3*d1*d1 + c2*c3*c1*d3*d3*d3 - c2*c3*c1*d3*d3*d1 - c2*c3*c1*d3*d1*d1 + c2*c3*c1*d1*d1*d1 - d2*d2*d3*c2*c3*c1
        - 6*d3*c1*c1*d2*d2 - 2*c3*c3*c3*c3*d1*d1*d2 - 2*d2*d2*c3*c3*c3*c3*d1
        - d2*d2*c3*c3*d1*c2*c2 - 4*d2*d2*c3*c3*d3 + d2*d2*c3*c3*c3*d1*c2*c1 - 2*d2*d2*d3*c2*c2
        + 4*d3*d3*d2*c2*c2 - d3*d3*d2*c2*c3*c1 + 6*d3*d3*d2*c1*c1 - 8*d3*d1*d2*c2*c2
        + 6*d3*d1*d2*c2*c3*c1 - 8*d3*d1*d2*c1*c1 + 4*d1*d1*d2*c2*c2 + 14*d2*d2*c3*c3*d1 + 48*d3*d1*d2 + 4*d3*d3*d1*c1*c1 - 2*d1*d1*d3*c1*c1 - 8*d1*d1*d1 + 8*d3*d3*d3;
    // x^0 coefficient
    coeffs[0] = ((d1+d2-d3)*(d1+d2-d3) - d2*c3*c3*d1 );
    coeffs[0] = coeffs[0]*coeffs[0];

    return coeffs;
}

// helper function for @ref fourPointPose computing the coefficients of the B matrix
static inline TooN::Vector<3> getBCoeffs( int i, int j, int k, int l, const TooN::Vector<5> & v4, const TooN::Vector<5> & v5 ){
    TooN::Vector<3> coeff;
    coeff[0] = v4[i]*v4[j] - v4[k]*v4[l];
    coeff[1] = v4[i]*v5[j] + v5[i]*v4[j] - (v4[k]*v5[l] + v5[k]*v4[l]);
    coeff[2] = v5[i]*v5[j] - v5[k]*v5[l];
    return coeff;
}

// contains the main computational part common to both of the high level functions
static bool fourPointSolver ( const std::vector<TooN::Vector<3> > & points, std::vector<TooN::Vector<3> > & myPixels, TooN::Vector<6> & distances, std::vector<TooN::Vector<2> > & length, const double angularError){
    TooN::Matrix<5> A;
    TooN::Vector<5> v4, v5;
    TooN::Matrix<7,3> B;
    TooN::Vector<3> bnull;
    TooN::Vector<5> t;
    TooN::Matrix<3> vecsCamera, vecsWorld;

    TooN::Vector<6> angles;

    // normalising scales for angle computation in next loop
    for(unsigned int i = 0; i < 4; i++){
        TooN::normalize(myPixels[i]);
    }

    int count = 0;
    for( unsigned int i = 0; i < 3; i++)
        for( unsigned int j = i+1; j < 4; j++, count++){
            TooN::Vector<3> diff = points[i] - points[j];
            distances[count] = diff * diff;
            angles[count] = 2* myPixels[i] * myPixels[j];
        }
	
	A = TooN::Zeros;
    A.slice<0,0,1,5>() = getACoeffs(angles[0], angles[1], angles[3], distances[0], distances[1], distances[3]).as_row();
    A.slice<1,0,1,5>() = getACoeffs(angles[0], angles[2], angles[4], distances[0], distances[2], distances[4]).as_row();
    A.slice<2,0,1,5>() = getACoeffs(angles[1], angles[2], angles[5], distances[1], distances[2], distances[5]).as_row();
    TooN::SVD<5> svdA(A);

    v4.as_row() = svdA.get_VT().slice<3,0,1,5>();
    v5.as_row() = svdA.get_VT().slice<4,0,1,5>();

    B.slice<0,0,1,3>() = getBCoeffs(4,2,3,3, v4, v5).as_row();
    B.slice<1,0,1,3>() = getBCoeffs(4,1,3,2, v4, v5).as_row();
    B.slice<2,0,1,3>() = getBCoeffs(4,0,3,1, v4, v5).as_row();
    B.slice<3,0,1,3>() = getBCoeffs(4,0,2,2, v4, v5).as_row();
    B.slice<4,0,1,3>() = getBCoeffs(3,1,2,2, v4, v5).as_row();
    B.slice<5,0,1,3>() = getBCoeffs(3,0,2,1, v4, v5).as_row();
    B.slice<6,0,1,3>() = getBCoeffs(2,0,1,1, v4, v5).as_row();

    TooN::SVD<7,3> svdB(B);
    bnull.as_row() = svdB.get_VT().slice<2,0,1,3>();

    double lambda, roh;
    if( bnull[1] != 0 ){          // lambda * roh != 0 => both are != 0
        double ratio = (bnull[0]/ bnull[1] + bnull[1] / bnull[2]) * 0.5;
        roh = 1.0 / (v4[0] * ratio + v5[0]);
        lambda = ratio * roh;
    } else if( bnull[2] != 0 ) { // roh != 0 => lambda == 0
        lambda = 0;
        roh = 1/v5[0];
    } else {                     // lambda != 0 and roh == 0
        roh = 0;
        lambda = 1/v4[0];
    }

    t = v4 * lambda + v5 * roh;
    double xx = (t[1]/t[0] + t[2]/t[1] + t[3]/t[2] + t[4]/t[3]) / 4;
    if( xx < 0)
        xx *= -1;

    bool valid = true;
    double x = sqrt( xx );
    length[0] = TooN::makeVector( x, -x); // possible distances to point 0
    valid &= computeDistances( xx, x, angles[0], distances[0], angularError, length[1]);
    valid &= computeDistances( xx, x, angles[1], distances[1], angularError, length[2]);
    valid &= computeDistances( xx, x, angles[2], distances[2], angularError, length[3]);
    return valid;
}

TooN::SE3<> fourPointPose( const std::vector<TooN::Vector<3> > & points, const std::vector<TooN::Vector<3> > & pixels, bool & valid, const double angularError ){
    double orientationTest = ((points[1] - points[0]) ^ (points[2] - points[0])) * (points[3] - points[0]);

    // normalising scales for angle computation in next loop
    std::vector<TooN::Vector<3> > myPixels(4);
    for(unsigned int i = 0; i < 4; i++)
        myPixels[i] = pixels[i];

    // resulting possible distances
    std::vector<TooN::Vector<2> > length(4);
    TooN::Vector<6> distances;
    valid = fourPointSolver( points, myPixels, distances, length, angularError);
    if( !valid )
        return TooN::SE3<>();

    // figure out the right lengths
    // brute force through all combinations, with some optimizations to stop early, if a better hypothesis exists
    std::vector<TooN::Vector<3> > testPoint(4);
    double minError = 1e10;
    int minIndex = -1;
    for(unsigned int h = 0; h < 16; h++){
        // first going through all possibilities where the point is in front of the camera
        // then switch to the other ones.
        testPoint[0] = myPixels[0] * length[0][!!(h&8)];
        testPoint[1] = myPixels[1] * length[1][!!(h&1)] * (h&8?-1:1);
        testPoint[2] = myPixels[2] * length[2][!!(h&2)] * (h&8?-1:1);
        testPoint[3] = myPixels[3] * length[3][!!(h&4)] * (h&8?-1:1);
        double error = 0;
        int count = 3; // skip the first 3 distances because they will not produce any error
        for( unsigned int i = 1; i < 3; i++)
            for( unsigned int j = i+1; j < 4; j++, count++){
                TooN::Vector<3> diff = testPoint[i] - testPoint[j];
                error += fabs(diff * diff - distances[count]);
            }
        if( minError > error ){
            double myOrientationTest = ((testPoint[1] - testPoint[0]) ^ (testPoint[2] - testPoint[0])) * (testPoint[3] - testPoint[0]);
            if(myOrientationTest * orientationTest > 0){
                minError = error;
                minIndex = h;
            }
        }
    }
    if(minIndex == -1){
        valid = false;
        return TooN::SE3<>();
    }

    // pixel directions extended to the right distance
    myPixels[0] *= length[0][!!(minIndex&8)];
    myPixels[1] *= length[1][!!(minIndex&1)] * (minIndex&8?-1:1);
    myPixels[2] *= length[2][!!(minIndex&2)] * (minIndex&8?-1:1);
    myPixels[3] *= length[3][!!(minIndex&4)] * (minIndex&8?-1:1);

    // absolute orientation
    return computeAbsoluteOrientation(points, myPixels);
}

// just a copy of the above code with modifications to the case search
TooN::SE3<> fourPointPoseFromCamera( const std::vector<TooN::Vector<3> > & points, const std::vector<TooN::Vector<3> > & pixels, bool & valid, const double angularError ){
    // normalising scales for angle computation in next loop
    std::vector<TooN::Vector<3> > myPixels(4);
    for(unsigned int i = 0; i < 4; i++)
        myPixels[i] = pixels[i];

    // resulting possible distances
    std::vector<TooN::Vector<2> > length(4);
    TooN::Vector<6> distances;
    valid = fourPointSolver( points, myPixels, distances, length, angularError);
    if( !valid )
        return TooN::SE3<>();

    // figure out the right lengths
    // brute force through all combinations, with some optimizations to stop early, if a better hypothesis exists
    std::vector<TooN::Vector<3> > testPoint(4);
    double minError = 1e10;
    int minIndex = -1;
    // we assume now that the point is in front of the camera and therefore the solution is the one with
    // positive x (which is the first in the vector)
    testPoint[0] = myPixels[0] * length[0][0];
    for(unsigned int h = 0; h < 8; h++){
        testPoint[1] = myPixels[1] * length[1][!!(h&1)];
        testPoint[2] = myPixels[2] * length[2][!!(h&2)];
        testPoint[3] = myPixels[3] * length[3][!!(h&4)];
        double error = 0;
        int count = 3; // skip the first 3 distances because they will not produce any error
        for( unsigned int i = 1; i < 3; i++)
            for( unsigned int j = i+1; j < 4; j++, count++){
                TooN::Vector<3> diff = testPoint[i] - testPoint[j];
                error += fabs(diff * diff - distances[count]);
            }
        if( minError > error ){
            minError = error;
            minIndex = h;
        }
    }
    if(minIndex == -1){
        valid = false;
        return TooN::SE3<>();
    }

    // pixel directions extended to the right distance
    myPixels[0] = testPoint[0];
    myPixels[1] *= length[1][!!(minIndex&1)];
    myPixels[2] *= length[2][!!(minIndex&2)];
    myPixels[3] *= length[3][!!(minIndex&4)];

    // absolute orientation
    return computeAbsoluteOrientation(points, myPixels);
}

}
