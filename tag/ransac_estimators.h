#ifndef TAG_RANSAC_ESTIMATORS_H
#define TAG_RANSAC_ESTIMATORS_H

#include <TooN/TooN.h>
#include <TooN/SVD.h>
#include <TooN/SymEigen.h>

namespace tag {

/// RANSAC estimator to compute an affine transformation between a set of 2D-2D correspondences.
/// The Correspondence datatype has to implement the following interface:
/// @code
/// struct Correspondence {
///   TooN::Vector<2> a;
///   TooN::Vector<2> b;
/// };
/// @endcode
/// The resulting transformation will map from a -> b.
/// @ingroup ransac
template <class Correspondence>
struct AffineHomography {
    /// the linear part of the resulting affine transformation
    TooN::Matrix<2> A;
    /// the translation part of the resulting affine transformation
    TooN::Vector<2> t;

    AffineHomography() { A[0][0] = A[0][1] = A[1][0] = A[1][1] = t[0] = t[1] = 0; }

    void estimate(const std::vector<Correspondence>& matches) {
        TooN::Matrix<> Ad(matches.size(), 3);
        TooN::Vector<> bx(matches.size()), by(matches.size());
        assert(matches.size() >= 3);
        for (unsigned int i=0; i<matches.size(); i++) {
            const Correspondence& m = matches[i];
            Ad[i].template slice<0,2>() =  m.a;
            Ad[i][2] = 1;
            bx[i] = m.b[0];
            by[i] = m.b[1];
        }
        //TooN::Matrix<3> ATA = Ad.T()*Ad;
        TooN::SVD<> svd(Ad);
        TooN::Vector<3> coefsX = svd.backsub(bx);
        TooN::Vector<3> coefsY = svd.backsub(by);
        A[0][0] = coefsX[0];
        A[0][1] = coefsX[1];
        t[0] = coefsX[2];
        A[1][0] = coefsY[0];
        A[1][1] = coefsY[1];
        t[1] = coefsY[2];
    }

    inline double getSqError(const Correspondence& m) const {
        TooN::Vector<2> disp = A*m.a + t - m.b;
        return disp*disp;
    }
};

/// RANSAC estimator to compute the fundamental matrix between two views based on 2D-2D correspondences.
/// The correspondences are coordinates in camera frame (after undistortion etc.)
/// The Correspondence datatype has to implement the following interface:
/// @code
/// struct Correspondence {
///   TooN::Vector<2> a;
///   TooN::Vector<2> b;
/// };
/// @endcode
/// @ingroup ransac
template <class Correspondence>
struct FundamentalMatrix {
    TooN::Matrix<3> F;
    TooN::Matrix<3> T;
    TooN::Vector<9> diag;

    void setScale(double scale, bool translate=true) {
        T[0][0] = T[1][1] = 1/scale;
        T[0][1] = T[1][0] = T[2][0] = T[2][1] = 0;
        T[0][2] = translate ? -0.5 : 0;
        T[1][2] = translate ? -0.5 : 0;
        T[2][2] = 1;
    }

    FundamentalMatrix() {}

    void estimate(const std::vector<Correspondence>& matches) {
        TooN::Matrix<> A(std::max(matches.size(),(unsigned)9),9);
        for (unsigned int i=0; i<matches.size(); i++) {
            TooN::Vector<3> a,b;
            a.template slice<0,2>() = matches[i].a;  a[2] = 1;
            b.template slice<0,2>() = matches[i].b;  b[2] = 1;
            a = T*a;
            b = T*b;
            A[i][0] = a[0]*b[0];
            A[i][1] = a[0]*b[1];
            A[i][2] = a[0];
            A[i][3] = a[1]*b[0];
            A[i][4] = a[1]*b[1];
            A[i][5] = a[1];
            A[i][6] = b[0];
            A[i][7] = b[1];
            A[i][8] = -1;
        }
        for (unsigned int i=matches.size(); i<9; i++) {
            for (int j=0; j<9; j++)
                A[i][j] = 0;
        }
        TooN::SVD<> svdA(A);
        diag = svdA.get_diagonal();
        TooN::Vector<9> f = svdA.get_VT()[8];
        for (int i=0; i<9; i++)
        F[i/3][i%3] = f[i];
        TooN::SVD<3> svdF(F);
        TooN::Vector<3> fdiag = svdF.get_diagonal();
        TooN::Matrix<3> DVT = svdF.get_VT();
        fdiag[2] = 0;
        for (int i=0; i<3; i++)
        DVT[i] *= fdiag[i]/fdiag[0];
        F = T.T() * svdF.get_U()*DVT * T;
    }

    inline double getSqError(const Correspondence& m) const {
        TooN::Vector<3> u;
        u[0] = m.b[0];
        u[1] = m.b[1];
        u[2] = 1;
        TooN::Vector<3> line = F*u;
        double dot = line[0]*m.a[0] + line[1]*m.a[1] + line[2];
        double a2 = line[0]*line[0];
        double b2 = line[1]*line[1];
        return (dot*dot)*1/(a2+b2);
    }
};

/// RANSAC estimator to compute ???
/// The Correspondence datatype has to implement the following interface:
/// @code
/// struct Correspondence {
///   TooN::Vector<2> a;
///   TooN::Vector<2> b;
/// };
/// @endcode
/// @ingroup ransac
template <class Correspondence>
struct DifferentialEssentialMatrix {
    TooN::Matrix<3> s;
    TooN::Vector<3> v,w;
    TooN::Vector<9> diag;

    DifferentialEssentialMatrix() {}

    void estimate(const std::vector<Correspondence>& matches) {
        TooN::Matrix<> A(std::max(matches.size(),(unsigned)9),9);
        for (unsigned int i=0; i<matches.size(); i++) {
            TooN::Vector<2> q=matches[i].a, u=matches[i].b-matches[i].a;
            A[i][0] = -u[1];
            A[i][1] = u[0];
            A[i][2] = u[1]*q[0]-u[0]*q[1];
            A[i][3] = q[0]*q[0];
            A[i][4] = 2*q[0]*q[1];
            A[i][5] = 2*q[0];
            A[i][6] = q[1]*q[1];
            A[i][7] = 2*q[1];
            A[i][8] = 1;
        }
        for (unsigned int i=matches.size(); i<9; i++) {
            for (int j=0; j<9; j++)
                A[i][j] = 0;
        }
        TooN::SVD<> svdA(A);
        diag = svdA.get_diagonal();
        TooN::Vector<9> e = svdA.get_VT()[8];
        // renormalize ?
        e *= 1.0/sqrt(e.template slice<0,3>() * e.template slice<0,3>());
        // Velocity
        v = e.template slice<0,3>();
        // Symmetric TooN::Matrix s
        s[0][0]=e[3];
        s[0][1]=s[1][0]=e[4];
        s[0][2]=s[2][0]=e[5];
        s[1][1]=e[6];
        s[1][2]=s[2][1]=e[7];
        s[2][2]=e[8];
        TooN::SymEigen<3> eigen(s);
        // We want e1 >= e2 >= e3, e1 >= 0, e3 <= 0
        TooN::Vector<3> lambda = eigen.get_evalues();
        lambda[0] = std::max(0.0,lambda[0]);
        lambda[2] = std::min(0.0,lambda[2]);
        // Project to special symmetric TooN::Matrix
        TooN::Vector<3> sigma = (TooN::make_Vector, 2*lambda[0] + lambda[1] - lambda[2],
                        lambda[0] + 2*lambda[1] + lambda[2],
                        -lambda[0] + lambda[1] + 2*lambda[2]);
        sigma /= 3;
        s = eigen.get_evectors().T() * diagmult(sigma,eigen.get_evectors());
        // recover w
        double lam = sigma[0] - sigma[2];
        double theta = acos(-sigma[1]/lam);
        TooN::Matrix<3> Ry1, Ry2, Rz;
        rotationY(theta/2 - M_PI/2, Ry1);
        rotationY(theta, Ry2);
        zero(Rz);
        Rz[0][1] = -1;
        Rz[1][0] = 1;
        Rz[2][2] = 1;
        TooN::Vector<3> diagLambda = (TooN::make_Vector,lam, lam,0);
        TooN::Vector<3> diag1 = (TooN::make_Vector,1,1,0);
        TooN::Matrix<3> V = eigen.get_evectors().T() * Ry1.T();
        TooN::Matrix<3> U = -1* V * Ry2;
        TooN::Vector<3> vhat[4] = { uncross(V*Rz*diagmult(diag1,V.T())),
                    uncross(V*Rz.T()*diagmult(diag1,V.T())),
                    uncross(U*Rz*diagmult(diag1,U.T())),
                    uncross(U*Rz.T()*diagmult(diag1,U.T())) };
        double dots[4];
        for (int k=0; k<4; k++)
        dots[k] = vhat[k]*v;
        int maxdot = std::max_element(dots, dots+4)-dots;
        switch (maxdot) {
        case 0: w = uncross(U*Rz*diagmult(diagLambda, U.T())); break;
        case 1: w = uncross(U*Rz.T()*diagmult(diagLambda, U.T())); break;
        case 2: w = uncross(V*Rz*diagmult(diagLambda, V.T())); break;
        case 3: w = uncross(V*Rz.T()*diagmult(diagLambda, V.T())); break;
        }
    }

    inline double getSqError(const Correspondence& m) const {
        TooN::Vector<3> q = (TooN::make_Vector, m.a[0], m.a[1], 1);
        TooN::Vector<3> u = (TooN::make_Vector, m.b[0]-m.a[0], m.b[1]-m.a[1], 0);
        TooN::Vector<3> top = v^q;
        TooN::Vector<3> bottom = s*q;
        return u*top + q*bottom;
    }

    TooN::Vector<3> recoverAngularVelocity() {
        TooN::SymEigen<3> eigen(s);
        TooN::Vector<3> sigma = eigen.get_evalues();
        double lambda = sigma[0] - sigma[2];
        double theta = acos(-sigma[1]/lambda);
        TooN::Matrix<3> Ry1, Ry2, Rz;
        rotationY(theta/2 - M_PI/2, Ry1);
        rotationY(theta, Ry2);
        zero(Rz);
        Rz[0][1] = -1;
        Rz[1][0] = 1;
        Rz[2][2] = 1;
        TooN::Vector<3> diagLambda = (TooN::make_Vector, lambda, lambda,0);
        TooN::Vector<3> diag1 = (TooN::make_Vector,1,1,0);
        TooN::Matrix<3> V = eigen.get_evectors().T() * Ry1.T();
        TooN::Matrix<3> U = -1* V * Ry2;
        TooN::Vector<3> vhat[4] = { uncross(V*Rz*diagmult(diag1,V.T())),
                    uncross(V*Rz.T()*diagmult(diag1,V.T())),
                    uncross(U*Rz*diagmult(diag1,U.T())),
                    uncross(U*Rz.T()*diagmult(diag1,U.T())) };
        double dots[4];
        for (int k=0; k<4; k++)
        dots[k] = vhat[k]*v;
        int maxdot = std::max_element(dots, dots+4)-dots;
        switch (maxdot) {
        case 0: w = uncross(U*Rz*diagmult(diagLambda, U.T())); break;
        case 1: w = uncross(U*Rz.T()*diagmult(diagLambda, U.T())); break;
        case 2: w = uncross(V*Rz*diagmult(diagLambda, V.T())); break;
        case 3: w = uncross(V*Rz.T()*diagmult(diagLambda, V.T())); break;
        }
        return w;
    }
};

} // namespace tag

#endif
