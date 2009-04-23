#include <tag/five_point.h>
#include <tag/helpers.h>

#include <iostream>
#include <fstream>

#include <TooN/TooN.h>

using namespace std;
using namespace TooN;

inline double rand_u() {
    return 2*((double) std::rand()/ RAND_MAX)-1;
}

template <int N> 
inline Matrix<N> normalize( Matrix<N> M ){
    M /= norm_inf(M);
    if(signbit(M(0,0)))
        M *= -1;
    return M;
}

typedef tr1::array<pair<Vector<3>, Vector<3> >,5 > InputData;

void test2Dfiledata(){
    tr1::array<pair<Vector<3>, Vector<3> >,5 > data;
    for(int i = 0; i < 5; ++i) {
        cin >> data[i].first.slice<0,2>().ref() >> data[i].second.slice<0,2>().ref();
       data[i].first[2] = 1;
       data[i].second[2] = 1;
    
        cout << data[i].first << "\t" << data[i].second << endl;
    }
    cout << endl;
    vector<Matrix<3> > result = tag::five_point(data);

    cout << "results\n" << endl;

    for(int i = 0; i < result.size(); ++i)
        cout << result[i]/norm_inf(result[i]) << endl;
}

void testTransform(){
    Vector<6> params;
    cin >> params;

    SE3<> T(params);
    Matrix<3> E = tag::getEssentialMatrix(T);

    cout << "initial E\n" << E/norm_inf(E) << endl;

    InputData data;

    data[0].first = makeVector(1, 0, 0);
    data[1].first = makeVector(0, 1, 0);
    data[2].first = makeVector(0, 0, 1);
    data[3].first = makeVector(1, 1, 1);
    data[4].first = makeVector(-1, 1, 1);

    for(int i = 0; i < 5; ++i) {
        data[i].second = T * data[i].first;
        cout << data[i].second * E * data[i].first << endl;
    }

    vector<Matrix<3> > result = tag::five_point(data);

    cout << "results\n" << endl;

    for(int i = 0; i < result.size(); ++i)
        cout << result[i]/norm_inf(result[i]) << endl;
}

void testRandom(){
    Vector<6> params;
    InputData data;

    srand(0); // same seed always

    double error = 0;

    const int ITERATIONS = 1000000;
    const double scale = 100;

    //ofstream diffile("diffs.txt");

    for(int i = 0; i < ITERATIONS; ++i){
        params = scale * makeVector(rand_u(), rand_u(), rand_u(), rand_u(), rand_u(), rand_u()); 
        SE3<> T(params);

        data[0].first = scale * makeVector(rand_u(), rand_u(), rand_u());
        data[1].first = scale * makeVector(rand_u(), rand_u(), rand_u());
        data[2].first = scale * makeVector(rand_u(), rand_u(), rand_u());
        data[3].first = scale * makeVector(rand_u(), rand_u(), rand_u());
        data[4].first = scale * makeVector(rand_u(), rand_u(), rand_u());

        for(int i = 0; i < 5; ++i) {
            data[i].second = T * data[i].first;
        }

        Matrix<3> E = normalize(tag::getEssentialMatrix(T));

        vector<Matrix<3> > result = tag::five_point(data);

        int index = -1;
        double diff = 1e100;
        for(int i = 0; i < result.size(); ++i){
            double d = norm_inf(E - normalize(result[i]));
            if( d < diff ){
                diff = d;
                index = i;
            }
        }
        if(index == -1){
            cout << "no solution !" << endl;
            cout << params << endl;
        } else {
            error += diff;
            //diffile << diff << "\n";
        }
    }
    cout << endl;
    cout << "average error " << error / ITERATIONS << endl;
}

int main(int argc, char ** argv){

    testRandom();

    return 0;
}
