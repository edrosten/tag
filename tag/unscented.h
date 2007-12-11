#ifndef TAG_UNSCENTED_H
#define TAG_UNSCENTED_H

#include <TooN/Cholesky.h>
#include <TooN/helpers.h>

namespace tag {

/// @defgroup unscented Unscented Transform
/// Implementations of the unscented transform for different sigma point selections. The @ref sphericalUnscentedTransform
/// will use only N + 2 sigma points (N dimension of input), while the @ref unscentedTransform uses 2N + 1 sigma points with
/// a simpler algorithm for sigma points. The @ref sphericalUnscentedTransform incurrs slight overhead in computing the sigma points,
/// but for expensive functions the difference in number evaluated points can win.

namespace Internal {

template <class F, int N, class V, class M> inline void outer_product_upper_half(const TooN::FixedVector<N,V>& v, const double w, TooN::FixedMatrix<N,N,M>& m)
{
    for (int i=0; i<N; ++i)
        for (int j=i; j<N; ++j)
            F::eval(m[i][j], w * v[i] * v[j]);
}

template <class F, int N, class V1, class V2, class M> inline void outer_product_upper_half(const TooN::FixedVector<N,V1>& v1, const TooN::FixedVector<N,V2>& v2, const double w, TooN::FixedMatrix<N,N,M>& m)
{
    for (int i=0; i<N; ++i)
        for (int j=i; j<N; ++j)
            F::eval(m[i][j], w * (v1[i] * v1[j] + v2[i] * v2[j]));
}

}

template <int N, int M, class F, class V1, class V2, class M1, class M2>
void unscentedTransformSqrt(const TooN::FixedVector<N, V1>& x, const TooN::FixedMatrix<N,N,M1>& L, const F& f, TooN::FixedVector<M,V2>& newx, TooN::FixedMatrix<M,M,M2>& newP)
{
    static const double w0 = 1/3.0;
    static const double w1 = (1-w0)/(2 * N);
    static const double spread = sqrt(N/(1-w0));
    const TooN::Vector<M> y0 = f(x);
    newx = w0 * y0;
    Internal::outer_product_upper_half<TooN::util::Assign>(y0, w0, newP);
    for (int i=0; i<N; i++) {
        const TooN::Vector<N> dxi = spread * L.T()[i];
        const TooN::Vector<M> yi1 = f(x + dxi);
        const TooN::Vector<M> yi2 = f(x - dxi);
        newx += w1 * (yi1 + yi2);
        Internal::outer_product_upper_half<TooN::util::PlusEquals>(yi1, yi2, w1, newP);
    }
    Internal::outer_product_upper_half<TooN::util::PlusEquals>(newx, -1, newP);
    TooN::Symmetrize(newP);
}

template <int N, int M, class F, class V1, class V2, class M1, class M2, class M3>
void unscentedTransformSqrt(const TooN::FixedVector<N, V1>& x, const TooN::FixedMatrix<N,N,M1>& L, const F& f, TooN::FixedVector<M,V2>& newx, TooN::FixedMatrix<M,M,M2>& newP, TooN::FixedMatrix<M,N,M3>& cov)
{
    static const double w0 = 1/3.0;
    static const double w1 = (1-w0)/(2 * N);
    static const double spread = sqrt(N/(1-w0));
    const TooN::Vector<M> y0 = f(x);
    newx = w0 * y0;
    Internal::outer_product_upper_half<TooN::util::Assign>(y0, w0, newP);
    TooN::Zero(cov);
    for (int i=0; i<N; i++) {
        const TooN::Vector<N> dxi = spread * L.T()[i];
        const TooN::Vector<M> yi1 = f(x + dxi);
        const TooN::Vector<M> yi2 = f(x - dxi);
        newx += w1 * (yi1 + yi2);
        Internal::outer_product_upper_half<TooN::util::PlusEquals>(yi1, yi2, w1, newP);
        TooN::add_product(w1 * (yi1-yi2).as_col(), dxi.as_row(), cov);
    }
    Internal::outer_product_upper_half<TooN::util::PlusEquals>(newx, -1, newP);
    TooN::Symmetrize(newP);
}

template<int N, int M, class F, class V1, class V2, class M1, class M2>
void sphericalUnscentedTransformSqrt(const TooN::FixedVector<N, V1> & x, const TooN::FixedMatrix<N,N,M1> & L, const F & f, TooN::FixedVector<M, V2> & newx, TooN::FixedMatrix<M, M, M2> & newP){
    static const double w0 = 1/3.0;
    static const double w1 = (1 - w0)/(N+1);
    static TooN::Matrix<N+1, N> xarg;
    static bool init = false;
    if(!init){
        init = true;
        Zero(xarg);
        for(int i = 0; i < N; ++i)
            xarg[0][i] = xarg[1][i] = -1.0/sqrt((i+1)*(i+2)*w1);
        xarg(0,0) *= -1.0;
        for(int i = 1; i < N; ++i){
            xarg[i+1][i] = -(i+1)*xarg[i][i];
            for(int j = i+1; j < N; ++j)
                xarg[i+1][j] = xarg[i][j];
        }
    }
    const TooN::Vector<M> y0 = f(x);
    newx = w0 * y0;
    Internal::outer_product_upper_half<TooN::util::Assign>(y0, w0, newP);
    for(int i = 0; i < N; i+=2){
        const TooN::Vector<M> res = f(x + L * xarg[i]);
        const TooN::Vector<M> res2 = f(x + L * xarg[i+1]);
        newx += w1 * (res+res2);
        Internal::outer_product_upper_half<TooN::util::PlusEquals>(res, res2, w1, newP);
    }
    if(N % 2 == 0){
        const TooN::Vector<M> res = f(x + L * xarg[N]);
        newx += w1 * res;
        Internal::outer_product_upper_half<TooN::util::PlusEquals>(res, w1, newP);
    }
    Internal::outer_product_upper_half<TooN::util::PlusEquals>(newx, -1, newP);
    TooN::Symmetrize(newP);
}

template<int N, int M, class F, class V1, class V2, class M1, class M2, class M3>
void sphericalUnscentedTransformSqrt(const TooN::FixedVector<N, V1> & x, const TooN::FixedMatrix<N,N,M1> & L, const F & f, TooN::FixedVector<M, V2> & newx, TooN::FixedMatrix<M, M, M2> & newP, TooN::FixedMatrix<M,N,M3>& cov ){
    static const double w0 = 1/3.0;
    static const double w1 = (1 - w0)/(N+1);
    static TooN::Matrix<N+1, N> xarg;
    static bool init = false;
    if(!init){
        init = true;
        Zero(xarg);
        for(int i = 0; i < N; ++i)
            xarg[0][i] = xarg[1][i] = -1.0/sqrt((i+1)*(i+2)*w1);
        xarg(0,0) *= -1.0;
        for(int i = 1; i < N; ++i){
            xarg[i+1][i] = -(i+1)*xarg[i][i];
            for(int j = i+1; j < N; ++j)
                xarg[i+1][j] = xarg[i][j];
        }
    }
    Zero(cov);
    const TooN::Vector<M> y0 = f(x);
    newx = w0 * y0;
    Internal::outer_product_upper_half<TooN::util::Assign>(y0, w0, newP);
    for(int i = 0; i < N+1; ++i){
        const TooN::Vector<N> d = L * xarg[i];
        const TooN::Vector<M> res = f(x + d);
        newx += w1 * res;
        Internal::outer_product_upper_half<TooN::util::PlusEquals>(res, w1, newP);
        TooN::add_product(w1 * (res - y0).as_col(), d.as_row(), cov);
    }
    Internal::outer_product_upper_half<TooN::util::PlusEquals>(newx, -1, newP);
    TooN::Symmetrize(newP);
}


/// computes the unscented transform of the distribution given by @a x and @a P through the function @a f using
/// 2*N + 1 sigma points distributed along the eigenvectors of the covariance matrix. Results are written into
/// @a newx and @a newP.
/// @param[in] x mean of the original distribution
/// @param[in] P covariance of the original distribution
/// @param[in] f function to apply
/// @param[out] newx mean of the result distribution
/// @param[out] newP covariance of the result distribution
/// @ingroup unscented
template <int N, int M, class F, class V1, class V2, class M1, class M2>
void unscentedTransform(const TooN::FixedVector<N, V1>& x, const TooN::FixedMatrix<N,N,M1>& P, const F& f, TooN::FixedVector<M,V2>& newx, TooN::FixedMatrix<M,M,M2>& newP)
{
    TooN::Matrix<N> L;
    TooN::Cholesky<N>::sqrt(P,L);
    unscentedTransformSqrt(x, L, f, newx, newP);
}

/// computes the unscented transform of the distribution given by @a x and @a P through the function @a f using
/// 2*N + 1 sigma points distributed along the eigenvectors of the covariance matrix. It also computes the covariance
/// between the original and result distribution. Results are written into @a newx, @a newP and @a cov.
/// @param[in] x mean of the original distribution
/// @param[in] P covariance of the original distribution
/// @param[in] f function to apply
/// @param[out] newx mean of the result distribution
/// @param[out] newP covariance of the result distribution
/// @param[out] cov covariance between input and output
/// @ingroup unscented
template <int N, int M, class F, class V1, class V2, class M1, class M2, class M3>
void unscentedTransform(const TooN::FixedVector<N, V1>& x, const TooN::FixedMatrix<N,N,M1>& P, const F& f, TooN::FixedVector<M,V2>& newx, TooN::FixedMatrix<M,M,M2>& newP, TooN::FixedMatrix<M,N,M3>& cov)
{
    TooN::Matrix<N> L;
    TooN::Cholesky<N>::sqrt(P,L);
    unscentedTransformSqrt(x, L, f, newx, newP, cov);
}

/// computes the unscented transform of the distribution given by @a x and @a P through the function @a f using
/// N + 2 sigma points distributed as a spherical simplex. Results are written into @a newx and @a newP.
/// @param[in] x mean of the original distribution
/// @param[in] P covariance of the original distribution
/// @param[in] f function to apply
/// @param[out] newx mean of the result distribution
/// @param[out] newP covariance of the result distribution
/// @ingroup unscented
template <int N, int M, class F, class V1, class V2, class M1, class M2>
void sphericalUnscentedTransform(const TooN::FixedVector<N, V1>& x, const TooN::FixedMatrix<N,N,M1>& P, const F& f, TooN::FixedVector<M,V2>& newx, TooN::FixedMatrix<M,M,M2>& newP){
    TooN::Matrix<N> L;
    TooN::Cholesky<N>::sqrt(P,L);
    sphericalUnscentedTransformSqrt(x, L, f, newx, newP);
}

/// computes the unscented transform of the distribution given by @a x and @a P through the function @a f using
/// N + 2 sigma points distributed as a spherical simplex. It also computes the covariance
/// between the original and result distribution. Results are written into @a newx, @a newP and @a cov.
/// @param[in] x mean of the original distribution
/// @param[in] P covariance of the original distribution
/// @param[in] f function to apply
/// @param[out] newx mean of the result distribution
/// @param[out] newP covariance of the result distribution
/// @param[out] cov covariance between input and output
/// @ingroup unscented
template <int N, int M, class F, class V1, class V2, class M1, class M2, class M3>
void sphericalUnscentedTransform(const TooN::FixedVector<N, V1>& x, const TooN::FixedMatrix<N,N,M1>& P, const F& f, TooN::FixedVector<M,V2>& newx, TooN::FixedMatrix<M,M,M2>& newP, TooN::FixedMatrix<M,N,M3>& cov){
    TooN::Matrix<N> L;
    TooN::Cholesky<N>::sqrt(P,L);
    sphericalUnscentedTransformSqrt(x, L, f, newx, newP, cov);
}

template <int N, int M, class Ax, class AP, class AR, class F>
void unscentedKalmanUpdate(TooN::FixedVector<N,Ax>& x, TooN::FixedMatrix<N,N,AP>& P, const F& f, const TooN::FixedMatrix<M,M,AR>& R)
{
    TooN::Vector<M> v;
    TooN::Matrix<M> Pyy;
    TooN::Matrix<M,N> Pyx;
    unscentedTransform(x, P, f, v, Pyy, Pyx);
    TooN::Cholesky<M> chol(Pyy + R);
    const TooN::Matrix<M> Sinv = chol.get_inverse();
    x += Pyx.T() * (Sinv * v);
    P -= transformCovariance(Pyx.T(), Sinv);
}

}

#endif
