#include <tag/fourpointpose.h>
#include <tag/absorient.h>

#include <TooN/SVD.h>
#include <TooN/helpers.h>

namespace tag {

// simple helper function to compute larger of a quadratic equation x^2 + px + q = 0
// also returns if the result is real or not (in which case valid is set to false)
static inline double quadraticRoots( double p, double q, bool & valid ){
    double det = p*p*0.25 - q;
    if( det < 0 ){
        valid = false;
        return 0;
    }
    return -p*0.5 + sqrt(det);
}

// helper function for @ref fourPointPose computing the coefficients of the A matrix
static TooN::Vector<5> getACoeffs(double c1, double c2, double c3, double d1, double d2, double d3){
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
static TooN::Vector<3> getBCoeffs( int i, int j, int k, int l, const TooN::Vector<5> & v4, const TooN::Vector<5> & v5 ){
    TooN::Vector<3> coeff;
    coeff[0] = v4[i]*v4[j] - v4[k]*v4[l];
    coeff[1] = v4[i]*v5[j] + v5[i]*v4[j] - (v4[k]*v5[l] + v5[k]*v4[l]);
    coeff[2] = v5[i]*v5[j] - v5[k]*v5[l];
    return coeff;
}

TooN::SE3 fourPointPose( const std::vector<TooN::Vector<3> > & points, const std::vector<TooN::Vector<3> > & pixels, bool & valid ){
    TooN::Matrix<5> A;
    TooN::Vector<5> v4, v5;
    TooN::Matrix<7,3> B;
    TooN::Vector<3> bnull;
    TooN::Vector<5> t;
    TooN::Vector<4> x;
    TooN::Matrix<3> vecsCamera, vecsWorld;

    TooN::Vector<6> distances;
    TooN::Vector<6> angles;

    // normalising scales for angle computation in next loop
    std::vector<TooN::Vector<3> > myPixels(4);
    for(unsigned int i = 0; i < 4; i++){
        myPixels[i] = pixels[i];
        TooN::normalize(myPixels[i]);
    }

    int count = 0;
    for( unsigned int i = 0; i < 3; i++)
        for( unsigned int j = i+1; j < 4; j++, count++){
            if(i == j)
                continue;
            TooN::Vector<3> diff = points[i] - points[j];
            distances[count] = diff * diff;
            angles[count] = 2* myPixels[i] * myPixels[j];
        }

    TooN::Zero(A);
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

    valid = true;
    x[0] = sqrt( xx );  // distance to point 0
    x[1] = quadraticRoots( -x[0]*angles[0], xx - distances[0], valid );
    x[2] = quadraticRoots( -x[0]*angles[1], xx - distances[1], valid );
    x[3] = quadraticRoots( -x[0]*angles[2], xx - distances[2], valid );

    if( !valid )
        return TooN::SE3();

    // pixel directions extended to the right distance
    for(unsigned int i = 0; i < 4; i++){
        myPixels[i] *= x[i];
    }

    // absolute orientation
    return computeAbsoluteOrientation(points, pixels);
}

}
