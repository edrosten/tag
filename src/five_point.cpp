#include <tag/five_point.h>
#include <tag/stdpp.h>
#include <tag/helpers.h>

#include <TooN/helpers.h>
#include <TooN/gauss_jordan.h>
#include <TooN/SVD.h>

#include <algorithm>
#include "polynomial.h"

using namespace TooN;
using namespace std;
using namespace std::tr1;


namespace tag {

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

	Vector<AN+BN-1> ret = Zeros;

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

Matrix<3, 3, double, Reference::RowMajor> as_matrix(Vector<9>& v)
{
	return Matrix<3, 3, double, Reference::RowMajor>(&v[0]);
}


// This function finds some vectors spanning the null space of m, assuming
// that none of the rows of m are linearly dependent.
//
//  m is a short, fat matrix.
//  The Gauss-Jordan decomposition is [I|A]
//
//  [I|A] * [A^T|-I]^T = [0]
//
// Since [I|A] spans the same space as m, then:
//
// m * [A^T|-I]^T = [0]
//
// Therefore [A^T|I] spans the null space of m.
template<int R, int C, class P> Matrix<C-R, C, P> dodgy_null(Matrix<R, C, P> m)
{
	gauss_jordan(m);

	Matrix<C-R, C> null;
	null.template slice<0,0,C-R, R>() = m.T().template slice<R,0,C-R, R>();
	null.template slice<0,R, C-R, C-R>() = -Identity;

	return null;
}


vector<Matrix<3> > five_point(const array<pair<Vector<3>, Vector<3> >, 5> & points)
{
	//Equations numbers are given with reference to:
	// "An efficient Solution to the Five-Point Relative Pose Problem",
	// D. Nister, IEEE Tran. on Pat. Anal. and Mach. Intel.,  26(6) 2004.
	// Franctional equations refer to equations between numbered ones.

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
	// E~ is in the null space of Q, and so the real E consists of a
	// linear sum of the remaining 4 null space vectors.
	// See Eqn 10.

	Matrix<5, 9> Q = Zeros;
	for(int i=0; i < 5; i++)
		Q[i] = stack_points(points[i].second, points[i].first);

	Matrix<4,9> null = dodgy_null(Q);

	// According to Eqn 10:
	Vector<9> X = null[0];
	Vector<9> Y = null[1];
	Vector<9> Z = null[2];
	Vector<9> W = null[3];

	Matrix<10,20> R;
	build_matrix(X, Y, Z, W, R);


	//Columns are:                               
	//                                            |    poly in x      poly in y        poly in 1
	//     LEFT HAND SIDE                         |  ______________   _____________   _________________
	//                                            |  '            '   '           '   '               '
	// x^3 y^3 x^2y xy^2 x^2z x^2 y^2z y^2 xyz xy |  z^2x   zx   x    z^2y   zy   y   z^3   z^2   z   1 

	gauss_jordan(R);

	// The left side is now the Identity matrix matching Eqn 10 1/2. 
	//
	// Due to the careful ordering of the coefficients, performing the
	// following subtractions (<k>, <l>, <m>) makes all parts on the left hand
	// side become zero, leaving only the right hand side (polynomials in x, y
	// and 1).

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

	//Compute the Determinant of B using expansion along the bottom row.
	//Eqn 24--28. 
	//
	// Det = b31 (b12 b23 -b13 b22) - b32 ( b11 b23 - b13 b21) + b33 ( b11 b22 - b12 b21)
	//     = b13  * p1 - b23  * p2 + b33 * p3
	Vector<8> p1 = poly_mul(b_12, b_23) - poly_mul(b_13, b_22);
	Vector<8> p2 = poly_mul(b_13, b_21) - poly_mul(b_11, b_23);
	Vector<7> p3 = poly_mul(b_11, b_22) - poly_mul(b_21, b_12);
	
	//The polynomial is...
	Vector<11> n = poly_mul(p1, b_31) + poly_mul(p2, b_32) + poly_mul(p3, b_33);

	vector<double> roots = find_roots(n);
	vector<Matrix<3> > Es;

	for(unsigned i=0; i <roots.size(); i++)
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

std::vector<TooN::SE3<> > se3_from_E( const TooN::Matrix<3> & E ){
	// follows the computation in 
	// Recovering Baseline and Orientation from 'Esssential' Matrix
	// BKP Horn, Jan 1990
	const Vector<3> & ea = E.T()[0];
	const Vector<3> & eb = E.T()[1];
	const Vector<3> & ec = E.T()[2];

	// cofactor matrix (19), used both to find largest vector product and later for rotation
	Matrix<3> cf;
	cf[0] = eb ^ ec;
	cf[1] = ec ^ ea;
	cf[2] = ea ^ eb;

	// find largest vector product
	const Vector<3> norms = makeVector(norm_sq(cf[0]), norm_sq(cf[1]), norm_sq(cf[2]));
	const int max_index = (int)(max_element(&norms[0], &norms[0]+3) - &norms[0]);

	// calculate direction vector at proper length (18)
	const Vector<3> t = unit(cf[max_index]) * sqrt(0.5 * (E[0] * E[0] + E[1] * E[1] + E[2] * E[2]));
	// scaling for rotation matrix later in (24)
	const double s = 1/(t*t);
	// the two rotation matrices using t and -t (24)
	const SO3<> Ra( s*(cf.T() - getCrossProductMatrix(t) * E) );
	const SO3<> Rb( s*(cf.T() - getCrossProductMatrix(-t) * E) );

	// put all the solutions together
	vector<TooN::SE3<> > SE3s;
	SE3s.push_back(SE3<>(Ra, t));
	SE3s.push_back(SE3<>(Ra, -t));
	SE3s.push_back(SE3<>(Rb, t));
	SE3s.push_back(SE3<>(Rb, -t));

	return SE3s;
}


TooN::SE3<> optimize_epipolar(const std::vector<std::pair<TooN::Vector<3>, TooN::Vector<3> > > & points, const TooN::SE3<> & initial){
	static const Vector<3> X = makeVector(1,0,0);

	// We represent E = [t]x R with t = Rn * [ 1 0 0 ]' and R = Rt
	// The parameterisation of the optimization is then
	// 3 parameters for Rt with left multiplication as update
	// 2 parameters for Rn with _right_ multiplication as update !!
	//   This allows as to have a minimal parameterization as 
	//   Rn * G_0 * [1 0 0]' = 0
	SO3<> Rt = initial.get_rotation();
	SO3<> Rn(X, initial.get_translation()); // rotation from X -> initial translation or Identity 

	int count = 0;
	WLS<5> wls;
	do {
		wls.clear();
		const Matrix<3> C = getCrossProductMatrix(Rn * X);
		for(unsigned i = 0; i < points.size(); ++i){
			Vector<5> J;
			const Vector<3> LEFT = points[i].second * C;
			const Vector<3> RIGHT = Rt * points[i].first;
			J[0] = LEFT * Rt.generator_field(0, RIGHT);
			J[1] = LEFT * Rt.generator_field(1, RIGHT);
			J[2] = LEFT * Rt.generator_field(2, RIGHT);
			J[3] = points[i].second * getCrossProductMatrix(Rn * Rn.generator_field(1, X)) * RIGHT;
			J[4] = points[i].second * getCrossProductMatrix(Rn * Rn.generator_field(2, X)) * RIGHT;
			const double e = 0 - LEFT * RIGHT;
			wls.add_mJ(e, J);
		}
		wls.compute();
		Rt = SO3<>::exp(wls.get_mu().slice<0,3>()) * Rt;
		Rn = Rn * SO3<>::exp(makeVector(0, wls.get_mu()[3], wls.get_mu()[4]));
		++count;
	} while(norm_sq(wls.get_mu()) > 1e-15 && count < 10);
	return SE3<>(Rt, Rn * X);
}


namespace {
double sq(double x)
{
	return x*x;
}
}
// Define a line l = [a b c] so that [x y 1] . [a b c] = 0
//
// The line equation can be written as:
//
// y = -ax/b -c/b
//
// And written in vector form:
//
// [ x ]   [0   ]              [ b ]
// [ y ] = [-c/b]  + x (1/b) * [ -a]
//
//The line normal therefore is v = [a b]
//
// r is the vector from [x0 y0] to any point on the line. The perpendicular distance to the
// line is |r . v| / |v|
//
// d = | (x - x0)*a + (y-y0) * b | / sqrt(a^2 + b^2)
//
// Rearranging and using ax + yb = -c:
//
// d = |ax0 + by0 + c| / sqrt(a*a + b*b)
double point_line_distance_squared(Vector<3> point, const Vector<3>& line)
{	
	//Normalize the point to [x0, y0, 1]
	point /= point[2];

	return sq(point * line) / (sq(line[0]) + sq(line[1]));
}

double point_line_distance(Vector<3> point, const Vector<3>& line)
{	
	//Normalize the point to [x0, y0, 1]
	point /= point[2];

	return point * line / sqrt(sq(line[0]) + sq(line[1]));
}


pair<double, double> essential_reprojection_errors_squared(const Matrix<3>& E, const Vector<3>& q, const Vector<3>& p)
{
	return make_pair(point_line_distance_squared(p, E*q), point_line_distance_squared(q, E.T()*p));
}

pair<double, double> essential_reprojection_errors(const Matrix<3>& E, const Vector<3>& q, const Vector<3>& p)
{
	return make_pair(point_line_distance(p, E*q), point_line_distance(q, E.T()*p));
}

}
