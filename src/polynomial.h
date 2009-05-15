#ifndef POLYNOMIAL_H
#define POLYNOMIAL_H
#include <vector>
#include <TooN/TooN.h>

std::vector<double> find_roots(const TooN::Vector<11>& v);


template<int N> double polyval(const TooN::Vector<N>& v, double x)
{
	//Polynomials are stored with the coefficient of zero in the first
	double val=0;
	for(int i=v.size()-1; i > 0; i--)
	{
		val += v[i];	
		val *= x;
	}

	val += v[0];
	return val;
}

#endif
