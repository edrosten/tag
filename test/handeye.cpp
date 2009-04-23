#include <tag/handeye.h>

#include <iostream>

using namespace std;
using namespace TooN;
using namespace tag;

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

int main(int argc, char ** argv){
    const int ITERATIONS = 100000;
    double total_error = 0;
    for(int t = 0; t < ITERATIONS; ++t){
        vector<SO3<> > a, b;
        
        SO3<> BC(rand_vect<3>(-2,2));
        SO3<> DA(rand_vect<3>(-2,2));
        
        for(int i = 0; i < 3; ++i ){
            SO3<> AB(rand_vect<3>(-2,2));
            a.push_back(AB);
            b.push_back(BC.inverse() * AB.inverse() * DA.inverse());
        }
        pair<SO3<>, SO3<> > xy = computeHandEye(a,b);
        double error = 0, e2 = 0, e3 = 0;
        for(unsigned i = 0; i < a.size(); ++i ){
            error +=  norm_sq((a[i] * xy.first * b[i] * xy.second).ln());
        }
        e2 = norm_sq((BC * xy.first.inverse()).ln());
        e3 = norm_sq((DA * xy.second.inverse()).ln());
        if( error > 1e-20 || e2 > 1e-20 || e3 > 1e-20 ){
            cout << "ups in rotation " << error << "\t" << e2 << "\t" << e3 << endl;
        }
        total_error += error;
    }
    cout << "SO3 average error\t" << total_error / ITERATIONS << endl;
    total_error = 0;
    for(int t = 0; t < ITERATIONS; ++t){
        vector<SE3<> > a, b;
        
        SE3<> BC(rand_vect<6>(-2,2));
        SE3<> DA(rand_vect<6>(-2,2));
        
        for(int i = 0; i < 3; ++i ){
            SE3<> AB(rand_vect<6>(-2,2));
            a.push_back(AB);
            b.push_back(BC.inverse() * AB.inverse() * DA.inverse());
        }
        pair<SE3<>, SE3<> > xy = computeHandEye(a,b);
        double error = 0;
        for(unsigned i = 0; i < a.size(); ++i ){
            error +=  norm_sq((a[i] * xy.first * b[i] * xy.second).ln());
        }
        if( error > 1e-20 ){
            cout << "ups in SE3 " << error << endl;
        }
        total_error += error;
    }
    cout << "SE3 average error\t" << total_error / ITERATIONS << endl;
    // for(double d = 1e-10; d < 1; d *= sqrt(10)){
        double d = 0.01;
        total_error = 0;
        for(int t = 0; t < ITERATIONS; ++t){
            vector<SE3<> > a, b;
            
            SE3<> BC(rand_vect<6>(-2,2));
            SE3<> DA(rand_vect<6>(-2,2));
            
            for(int i = 0; i < 4; ++i ){
                SE3<> AB(rand_vect<6>(-2,2));
                a.push_back(AB * SE3<>::exp(rand_vect<6>(-d,d)));
                b.push_back(BC.inverse() * AB.inverse() * DA.inverse() * SE3<>::exp(rand_vect<6>(-d,d)));
            }
            pair<SE3<>, SE3<> > xy = computeHandEye(a,b);
            double error = 0;
            for(unsigned i = 0; i < a.size(); ++i ){
                error +=  norm_sq((a[i] * xy.first * b[i] * xy.second).ln());
            }
            total_error += error;
        }
        cout << "pertupated SE3 with " << d << "\t" << total_error / ITERATIONS << endl;
    //}
    return 0;
}
