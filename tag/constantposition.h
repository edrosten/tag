#ifndef TAG_CONSTANTPOSITION_H
#define TAG_CONSTANTPOSITION_H

#include <TooN/se3.h>

namespace tag {

namespace ConstantPosition {

/// @defgroup constantpositiongroup Constant Position
/// contains State, Model and Measurement classes for a constant position model for use
/// with the KalmanFilter class.
/// @ingroup kalmanfiltergroup

/// The State class containing an SE3 to represent the filter state and the covariance
/// matrix of it.
/// @ingroup constantpositiongroup
class State {
public:
    State(void){
        reset();
    }

    void reset(void){
        pose = TooN::SE3<>();
        covariance = TooN::Identity;
    }

    static const int STATE_DIMENSION = 6;
    TooN::SE3<> pose;
    TooN::Matrix<STATE_DIMENSION> covariance;
};

/// operator to print out instances of State in a useable manner.
/// @ingroup constantpositiongroup
template <class O>
O & operator<< (O & os , const State & st){
    os << st.pose.ln() << st.pose.inverse().get_translation();
    return os;
}

/// The Model class implementing a constant position model. It will
/// only update the covariance based on time passing, and update
/// the state correctly from a measurement.
/// @ingroup constantpositiongroup
class Model {
public:
    /// describes the process noise as independent for each dimension of the state
    TooN::Vector<State::STATE_DIMENSION> sigma;
    /// the jacobian of the process modell, here the identity
    TooN::Matrix<State::STATE_DIMENSION> jacobian;
    /// the actual process noise matrix returned from the associated funcion
    TooN::Matrix<State::STATE_DIMENSION> noise;

    Model(void){
        sigma = TooN::Zeros;
        noise = TooN::Zeros;
        jacobian = TooN::Identity;
    }

    /// Jacobian has pos, rot in this order
    TooN::Matrix<State::STATE_DIMENSION> & getJacobian(const State & state, double dt){
        return jacobian;
    }

    void updateState( State & state, const double dt ){
    }

    TooN::Matrix<State::STATE_DIMENSION> & getNoiseCovariance( double dt ){
        for(unsigned int i = 0; i < 6; i++){
            noise(i,i) = dt * sigma[i];
        }
        return noise;
    }

    void updateFromMeasurement( State & state, const TooN::Vector<State::STATE_DIMENSION> & innovation ){
        state.pose = TooN::SE3<>::exp(innovation) * state.pose;
    }
};

} // namespace ConstantPosition

} // namespace tag

#endif
