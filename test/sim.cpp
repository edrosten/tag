#include <iostream>
#include <vector>

#include <tag/absorient.h>

using namespace std;
using namespace TooN;
using namespace tag;

typedef vector<Vector<3> > Points;

double inline rand_u(){
    return (double)rand()/RAND_MAX;
}

double inline rand_u(double from, double to){
    return from + (to - from)*rand_u();
}

template <int N>
Vector<N> rand_vect(double from, double to){
    Vector<N> v;
    for(int i = 0; i < N; ++i)
        v[i] = rand_u(from, to);
    return v;
}

int main( int arg, char ** argv ){

    const int N = arg > 1 ? atoi(argv[1]) : 10;

    Points a(N), b(N), c(N);
    
    SIM3<> pose(makeVector(1,2,3,2,1,0.5, 0.37));
    
    for( unsigned i = 0; i < N; ++i ){
        a[i] = rand_vect<3>(-10, 10);
        b[i] = pose * a[i];
        c[i] = b[i] + rand_vect<3>(-0.01, 0.01);
    }

    SIM3<> sb = computeSimilarity(a,b);
    SIM3<> sc = computeSimilarity(a,c);

    double eb = 0, ec = 0;
    for( unsigned i = 0; i < N; ++i){
        eb += norm_sq(b[i] - sb * a[i]);
        ec += norm_sq(c[i] - sc * a[i]);
    }
    eb = sqrt(eb/N);
    ec = sqrt(ec/N);
    
    cout << "test\t" << pose.ln() << "\n";
    cout << "exact\t" << sb.ln()  << "\t" << eb << "\n"
         << "noisy\t" << sc.ln() << "\t" << ec << endl;

    return 0;
}
