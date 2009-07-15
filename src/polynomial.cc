#include <TooN/TooN.h>
#include <utility>
#include <vector>
#ifndef WIN32
#include <tr1/tuple>
#include <tr1/array>
#else
#include <tuple>
#include <array>
#endif

#include "polynomial.h"

using namespace TooN;
using namespace std;
using namespace std::tr1;

#ifdef WIN32
inline bool signbit( const double & d ){
	return d < 0.0;
}
#endif

template<int AN, int BN> Vector<AN+BN-1> poly_mul(const Vector<AN>& a, const Vector<BN>& b)
{
	//Polynomials are stored with the coefficient of zero in the first

	Vector<AN+BN-1> ret = Zeros;

	for(int ai=0; ai < a.size(); ai++)
		for(int bi=0; bi <b.size(); bi++)
			ret[ai+bi] += a[ai] * b[bi];
			
	return ret;
}
template<int N> Vector<N-1> polydiff(const Vector<N>& v)
{
	Vector<N-1> ret;

	for(int i=1; i < N; i++)
		ret[i-1] = v[i] * i;

	return ret;
}

template<int NumNum, int NumDenom>
pair<Vector<NumNum-NumDenom+1>, Vector<NumDenom-1> > poly_div(Vector<NumNum> num, const Vector<NumDenom>& denom)
{
	static const int NumQuot = NumNum-NumDenom+1;
	static const int NumRem = NumDenom - 1;

	Vector<NumQuot> quotient;

	for(int i=NumQuot-1; i >=0 ; i--)
	{
		double scale = num[i+denom.size()-1] / denom[denom.size()-1];
		quotient[i] = scale;

		for(int j=0; j < denom.size(); j++)
			num[i+j] -= scale * denom[j];
	}
	
	Vector<NumRem> remainder = num.template slice<0,NumRem>();
	return make_pair(quotient, remainder);	
}

template<int N, int M>
tuple<Vector<N>, Vector<M> > neg_second(const pair<Vector<N>, Vector<M> >& r)
{
	return make_pair(r.first, -r.second);
}

struct SturmChain10
{
	Vector<11> f10;
	Vector<10> f9;
	Vector<9 > f8;
	Vector<8 > f7;
	Vector<7 > f6;
	Vector<6 > f5;
	Vector<5 > f4;
	Vector<4 > f3;
	Vector<3 > f2;
	Vector<2 > f1;
	Vector<1 > f0;

	array<Vector<2>, 11> q;

	SturmChain10(const Vector<11>& p)
	{
		f10 = p;
		f9  = polydiff(p);

		tie(q[10], f8) = neg_second(poly_div(f10, f9));
		tie(q[ 9], f7) = neg_second(poly_div(f9,  f8));
		tie(q[ 8], f6) = neg_second(poly_div(f8,  f7));
		tie(q[ 7], f5) = neg_second(poly_div(f7,  f6));
		tie(q[ 6], f4) = neg_second(poly_div(f6,  f5));
		tie(q[ 5], f3) = neg_second(poly_div(f5,  f4));
		tie(q[ 4], f2) = neg_second(poly_div(f4,  f3));
		tie(q[ 3], f1) = neg_second(poly_div(f3,  f2));
		tie(q[ 2], f0) = neg_second(poly_div(f2,  f1));
	}

	// The sturm chain is constructed so that:
	//
	//  f_n   = P
	//  f_n-1 = P'
	//  f_n-2 = -REM(f_n, f_n-1)
	//  ...
	//
	// Simply from division:
	// 
	//  f_i = QUOT(f_i, f_i-1) * f_i-1 + REM(f_i, r_i-1)
	//
	//  From the sequence, then:
	//
	//  f_i = Q_i * f_i-1 - f_i-2
	//
	// Q_i is always of degree 1, so the number of sign changes
	// can be evaluated using the recursion above efficiently.
	tuple<int, double> changes(double x) const
	{
		int changes=0;
		double v_prev_2 = polyval(f0, x);
		double v_prev_1 = polyval(f1, x);

		if(signbit(v_prev_1) != signbit(v_prev_2))
			changes++;

		for(int i=2; i < 11; i++)
		{
			double v = polyval(q[i], x) * v_prev_1 - v_prev_2;

			if(signbit(v) != signbit(v_prev_1))
				changes++;
			v_prev_2 = v_prev_1;
			v_prev_1 = v;
		}
		return make_pair(changes, v_prev_1);
	}
	
	// At infinity, all aditive terms disappear
	int changes_at_infinity() const
	{
		int signs = 0;

		bool s_0 = signbit(f0[0]);
		bool s_prev = signbit(f1[1]);

		if(s_0 != s_prev)
			signs++;

		for(int i=2; i < 11; i++)
		{
			bool s = s_prev != (bool)signbit(q[i][1]);

			if(s != s_prev)
				signs++;
			
			s_prev = s;
		}
		return signs;	
	}

	int changes_at_neg_infinity() const
	{
		int signs = 0;

		bool s_0 = signbit(f0[0]);
		bool s_prev = !signbit(f1[1]);

		if(s_0 != s_prev)
			signs++;

		for(int i=2; i < 11; i++)
		{
			bool s = s_prev != (bool)!signbit(q[i][1]);

			if(s != s_prev)
				signs++;
			
			s_prev = s;
		}
		return signs;	
	}

	double operator()(double x) const
	{
		return polyval(f10, x);
	}

	double deriv(double x) const
	{
		return polyval(f9, x);
	}
};

template<class F> double polish_root_bisection(double lower, double upper, double lower_val, double upper_val, const F& f)
{
	const double zeps = 1e-20;
	const double eps  = 1e-10;
	static const int maxits = 50;

	int it;
	
	for(it = 0; fabs(upper-lower) / (fabs(lower) + fabs(upper)) > eps && fabs(upper-lower) > zeps && it < maxits; it++)
	{
		double midpoint = (upper + lower)/2;
		double mid_val = f(midpoint);

		if((bool)signbit(mid_val) == (bool)signbit(upper_val))
		{
			upper_val = mid_val;
			upper = midpoint;
		}
		else
		{
			lower_val = mid_val;
			lower = midpoint;
		}
	}

	cerr << "Converged in " << it << " iterations\n";
	
	return (upper + lower)/2;
}

template<class F> tuple<double, double> polish_root_newton(double lower, double upper, double lower_val, double upper_val, const F& f)
{
	const double zeps = 1e-20;
	const double eps  = 1e-15;
	static const int maxits = 10;

	double r =  (lower + upper)/2;		//Initial root guess
	double fr  = f(r);
	double fdr = f.deriv(r);

	double dxold = upper-lower;
	double dx = dxold;
	
	int it;
	
	for(it=0; it < maxits; it++)
	{
		double newton_step = -fr / fdr;
		
		//Check suitability of newton step.
		//Must be within bounds and reasonably quick convergence
		if(r + newton_step > upper || r + newton_step < lower || 2*fabs(newton_step) > fabs(dxold))
		{
			//Unsuitable, so take a bisection step
			dxold = dx;
			dx = (upper - lower) / 2;
			r  = lower + dx;  //r = (xl+xh)/2
		}
		else
		{
			r += newton_step;

			dxold = dx;
			dx = newton_step;
		}

		//Check for convergence
		if( fabs(dx / r) < eps || fabs(dx) < zeps)
			return make_pair(r, fr);

		//Update
		fr = f(r);
		fdr = f.deriv(r);

		if( (bool)signbit(fr) == (bool)signbit(upper_val))
		{
			upper_val = fr;
			upper = r;
		}
		else
		{
			lower_val = fr;
			lower = r;
		}
	}
	
	return make_pair(r, fr);
}




vector<double> find_roots(const Vector<11>& v)
{
	//This is a coarse root-finding routine. No effort is made to 
	//track down really reluctant roots, they are simply discarded.
	//Easy roots are found to a decent accuracy.

	const double zeps = 1e-20;  
	const double eps  = 1e-10;  
	const double feps = 1e-4;	//Epsilon for the function value
	
	//Roots above this are impossible to determine accurately
	//since the polynomial is so steep that even machine precision
	//errors in the root can lead to vast errors.
	const double highest_root = 100;

	SturmChain10 s(v);

	vector<double> roots;
	
	int at_n_inf = s.changes_at_neg_infinity();
	int at_p_inf = s.changes_at_infinity();
	
	int num_roots = at_n_inf - at_p_inf;

	//Try to find an upper bound for the highest root
	double upper=1, upper_val;

	for(;;)
	{
		int c;
		tie(c, upper_val) = s.changes(upper);

		
		if(at_n_inf - c < num_roots && upper < highest_root)
			upper*=2;
		else
			break;
	}
	
	//Try to find a lower bound for the lowest roots
	double lower=-1, lower_val;
	for(;;)
	{
		int c;
		tie(c, lower_val) = s.changes(lower);
		
		if(c - at_p_inf < num_roots && lower > -highest_root)
			lower*=2;
		else
			break;
	}


	//Now we have bounds for all roots and some assosciated polynomial
	//values.


	vector<tuple<double, double, int, int, double, double> > intervals;
	intervals.push_back(make_tuple(lower, upper, at_n_inf, at_p_inf, lower_val, upper_val));

	while(!intervals.empty())
	{
		int at_lower, at_upper;

		tie(lower, upper, at_lower, at_upper, lower_val, upper_val) = intervals.back();
		intervals.pop_back();

		double midpoint = (upper + lower)/2;
		
		//If the bracket is good, then continue, otherwise discard it.
		if(fabs((upper-lower) / midpoint) > eps && fabs(upper - lower) > zeps)
		{
			int at_midpoint;
			double midpoint_val;

			tie(at_midpoint, midpoint_val) = s.changes(midpoint);


			if(at_midpoint == at_lower)
				intervals.push_back(make_tuple(midpoint, upper, at_midpoint, at_upper, midpoint_val, upper_val));
			else if(at_midpoint == at_upper)
				intervals.push_back(make_tuple(lower, midpoint, at_lower, at_midpoint, lower_val, midpoint_val));
			else
			{
				//We've split the interval

				//Check to see if the lower interval has bracketed a single 
				//root
				if(at_lower - at_midpoint == 1)
				{
					double root, rootval;

					tie(root, rootval) = polish_root_newton(lower, midpoint, lower_val, midpoint_val, s);
					
					//Discard poor-quality roots
					if(fabs(rootval) < feps)
						roots.push_back(root);
				}
				else
					intervals.push_back(make_tuple(lower, midpoint, at_lower, at_midpoint, lower_val, midpoint_val));

				//Do the same with the upper interval
				if(at_midpoint - at_upper == 1)
				{
					double root, rootval;

					tie(root, rootval) = polish_root_newton(midpoint, upper, midpoint_val, upper_val, s);

					//Discard poor-quality roots
					if(fabs(rootval) < feps)
						roots.push_back(root);
				}
				else
					intervals.push_back(make_tuple(midpoint, upper, at_midpoint, at_upper, midpoint_val, upper_val));
			}
		}
	}

	return roots;
}


