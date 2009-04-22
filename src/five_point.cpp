#include <tag/five_point.h>
#include <tag/stdpp.h>

#include <TooN/helpers.h>
#include <TooN/gauss_jordan.h>
#include <TooN/SVD.h>

using namespace TooN;
using namespace std;
using namespace std::tr1;

extern "C" {
#include "solve.h"
}

namespace tag {

vector<double> get_roots(const Vector<11> & p){
	// uses Graphics Gem I code 
	// for details see here:
	// http://tog.acm.org/GraphicsGems/

	poly	sseq[10];
	Vector<11,double,Reference> v(sseq[0].coef);
	v = p;

	int num_poly = buildsturm(10, sseq);
	int atmin, atmax;
	int nroots = numroots(num_poly, sseq, &atmin, &atmax);

	vector<double> roots(nroots);
	if(nroots == 0)
		return roots;

	double min = -1.0;
	int nchanges = numchanges(num_poly, sseq, min);
	for (int i = 0; nchanges != atmin && i != MAXPOW; i++) { 
		min *= 10.0;
		nchanges = numchanges(num_poly, sseq, min);
	}

	if (nchanges != atmin) {
			cerr << "solve: unable to bracket all negative roots" << endl;
			atmin = nchanges;
	}

	double max = 1.0;
	nchanges = numchanges(num_poly, sseq, max);
	for (int i = 0; nchanges != atmax && i != MAXPOW; i++) { 
		max *= 10.0;
		nchanges = numchanges(num_poly, sseq, max);
	}

	if (nchanges != atmax) {
			cerr << "solve: unable to bracket all positive roots" << endl;
			atmax = nchanges;
	}

	nroots = atmin - atmax;

	/*
	 * perform the bisection.
	 */
	sbisect(num_poly, sseq, min, max, atmin, atmax, &roots[0]);
}

void build_matrix(const Vector<9>& X, const Vector<9>& Y, const Vector<9>& Z, const Vector<9>& W, Matrix<10,20>& R);

Vector<9> stack_points(const Vector<3>&p, const Vector<3>& q)
{
	return makeVector(
		q[0] *p[0],
		q[1] *p[0],
		q[2] *p[0],

		q[0] *p[1],
		q[1] *p[1],
		q[2] *p[1],

		q[0] *p[2],
		q[1] *p[2],
		q[2] *p[2]
	);
}

template<int AN, int BN> Vector<AN+BN-1> poly_mul(const Vector<AN>& a, const Vector<BN>& b)
{
	//Polynomials are stored with the coefficient of zero in the first

	Vector<AN+BN-1> ret = Zero;

	for(int ai=0; ai < a.size(); ai++)
		for(int bi=0; bi <b.size(); bi++)
			ret[ai+bi] += a[ai] * b[bi];
			
	return ret;
}

template<int N> double polyval(const Vector<N>& v, double x)
{
	double val = 0;
	for(int i=v.size()-1; i > 0; i--)
	{
		val += v[i];	
		val *= x;
	}

	val += v[0];
	return val;
}

template<int N> Vector<N-1> poly_diff(const Vector<N> & v)
{
	Vector<N-1> ret;
	for(int i = 1; i < N; ++i){
		ret[i-1] = v[i]*i;
	}
	return ret;
}

template<int N> pair<Vector<2>, Vector<N-1> > poly_div(const Vector<N> & num, const Vector<N-1> & denom)
{
	Vector<2> f;
	Vector<N-1> r;
	
	
}

Matrix<3, 3, double, Reference::RowMajor> as_matrix(Vector<9>& v)
{
	return Matrix<3, 3, double, Reference::RowMajor>(&v[0]);
}

vector<Matrix<3> > five_point(array<pair<Vector<3>, Vector<3> >, 5> points)
{
	//Equations numbers are given with reference to:
	// "An efficient Solution to the Five-Point Relative Pose Problem",
	// D. Nister, IEEE Tran. on Pat. Anal. and Mach. Intel.,  26(6) 2004.

	//Given a pair of corresponding points, p, q, the
	//epipolar constraint gives:
	//
	// p E q = 0
	//
	// This can be arranged in to:
	//
	// q~ E~ = 0
	//
	// Where E~ = [E11 E12 E13 .... ]  (ie E rastered in row-major order)
	// and q~ = [ q1p1 q2p1 q3p1 q1p2 q2p2 q3p2 q1p3 q2p3q3p3 ]
	// See Eqn 4, 7, 8, 9

	// Since there are 5 points, All 5 q~'s can be stacked giving Q, and
	//
	// Q E~ = [0 0 0 0 0]'
	//
	// E~ is in the null space of Q, and so the real E# consists of a
	// linear sum of the remaining 4 null space vectors.
	// See Eqn 10.

	Matrix<9, 9> Q = Zero;
	for(int i=0; i < 5; i++)
		Q[i] = stack_points(points[i].first, points[i].second);

	SVD<9, 9> svd_Q(Q);

	//The null-space it the last 4 rows of svd_Q.get_VT()
	//
	// According to Eqn 10:
	Vector<9> X = svd_Q.get_VT()[5];
	Vector<9> Y = svd_Q.get_VT()[6];
	Vector<9> Z = svd_Q.get_VT()[7];
	Vector<9> W = svd_Q.get_VT()[8];

	
	Matrix<10,20> R;
	build_matrix(X, Y, Z, W, R);


	//Columns are:
	gauss_jordan(R);


	//Build the B matrix (Eqn 13 1/2)
	//Polynomials of degree N are stored in a Vector<N+1> with Vector[0] as the coefficient of 1.
	//
	//These come from the right hand part of R
	// 10     11   12   13     14  15   16    17   18   19
	// z^2x   zx   x    z^2y   zy   y   z^3   z^2   z   1 

	//<k> = <e> - z<f>   e == Row 4, f == Row 5
	//<l> = <g> - z<h>   g == Row 6, h == Row 7
	//<m> = <i> - z<j>   i == Row 8, j == Row 9

	/*
        [ b11 b12 b13 ]
	B = [ b21 b22 b23 ] 
        [ b31 b32 b33 ]
	*/
	
	//b_11 is <k> coefficients of x
	//b_21 is <l> coefficients of x
	//b_31 is <m> coefficients of x
	Vector<4> b_11 = makeVector(R[4][12], R[4][11], R[4][10], 0) - makeVector(0, R[5][12], R[5][11], R[5][10]);
	Vector<4> b_21 = makeVector(R[6][12], R[6][11], R[6][10], 0) - makeVector(0, R[7][12], R[7][11], R[7][10]);
	Vector<4> b_31 = makeVector(R[8][12], R[8][11], R[8][10], 0) - makeVector(0, R[9][12], R[9][11], R[9][10]);
	
	//b_12 is <k> coefficients of y
	//b_22 is <l> coefficients of y
	//b_32 is <R> coefficients of y
	Vector<4> b_12 = makeVector(R[4][15], R[4][14], R[4][13], 0) - makeVector(0, R[5][15], R[5][14], R[5][13]);
	Vector<4> b_22 = makeVector(R[6][15], R[6][14], R[6][13], 0) - makeVector(0, R[7][15], R[7][14], R[7][13]);
	Vector<4> b_32 = makeVector(R[8][15], R[8][14], R[8][13], 0) - makeVector(0, R[9][15], R[9][14], R[9][13]);

	//b_12 is <k> coefficients of 1
	//b_22 is <l> coefficients of 1
	//b_32 is <R> coefficients of 1
	Vector<5> b_13 = makeVector(R[4][19], R[4][18], R[4][17], R[4][16], 0) - makeVector(0, R[5][19], R[5][18], R[5][17], R[5][16]);
	Vector<5> b_23 = makeVector(R[6][19], R[6][18], R[6][17], R[6][16], 0) - makeVector(0, R[7][19], R[7][18], R[7][17], R[7][16]);
	Vector<5> b_33 = makeVector(R[8][19], R[8][18], R[8][17], R[8][16], 0) - makeVector(0, R[9][19], R[9][18], R[9][17], R[9][16]);

	//Compute the Determinant of B using expansion by minors
	//Not quite Eqn 24--28. Nister uses col/row not row/col
	//
	// Det = b13 (b21 b32 -b31 b22) - b23 ( b11 b32 - b31 b12) + b33 ( b11 b22 - b21 b12)
	//     = b13  * p1 - b23  * p2 + b33 * p3
	Vector<8> p1 = poly_mul(b_12, b_23) - poly_mul(b_13, b_22);
	Vector<8> p2 = poly_mul(b_13, b_21) - poly_mul(b_11, b_23);
	Vector<7> p3 = poly_mul(b_11, b_22) - poly_mul(b_21, b_12);
	
	//The polynomial is...
	Vector<11> n = poly_mul(p1, b_31) + poly_mul(p2, b_32) + poly_mul(p3, b_33);

	vector<double> roots = get_roots(n);
	vector<Matrix<3> > Es;

	for(int i=0; i <roots.size(); i++)
	{
		double z = roots[i];
		
		//Solve the linear system in x, y for the forst two rows of b
		// [b11 b12 b13] [x]   [ 0 ]
		// [b21 b22 b23] [y] = [ 0 ]
		//
		//  Rearrange to give:
		//  
		// [b11 b12] [x]    [ b13 ]
		// [b21 b22] [y] = -[ b23 ]
		//
		// Solve and it gives:
		// Eqn 28. 
		double x = polyval(p1, z)/polyval(p3,z);
		double y = polyval(p2, z)/polyval(p3,z);

		Es.push_back(x * as_matrix(X) + y*as_matrix(Y) + z*as_matrix(Z) + as_matrix(W));
	}

	return Es;
}

}
