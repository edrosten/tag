#ifndef TAG_CONSTANTVELOCITY_H
#define TAG_CONSTANTVELOCITY_H

#include <TooN/se3.h>

namespace tag {

namespace ConstantVelocity {

/// @defgroup constantvelocitygroup Constant Velocity
/// contains State, Model and Measurement classes for a constant velocity model for use
/// with the KalmanFilter class.
/// @ingroup kalmanfiltergroup

/// The State class containing an SE3 for pose a Vector<3> for velocity and Vector<3> for angularVelocity
/// to represent the filter state and the covariance matrix of it.
/// @ingroup constantvelocitygroup
class State {
public:
    State(void){
        reset();
    }

    inline void reset(void){
        pose = TooN::SE3();
        TooN::Zero(angularVelocity);
        TooN::Zero(velocity);
        TooN::Identity(covariance);
    }

    inline void resetVelocity(void){
        TooN::Zero(angularVelocity);
        TooN::Zero(velocity);
    }

    static const int STATE_DIMENSION = 12;
    TooN::SE3 pose;
    TooN::Vector<3> angularVelocity;
    TooN::Vector<3> velocity;
    TooN::Matrix<STATE_DIMENSION> covariance;
};

/// operator to print out instances of State in a useable manner.
/// @ingroup constantvelocitygroup
template <class O>
inline O & operator<< (O & os , const State & st){
    os << st.pose.ln() << st.velocity << st.angularVelocity << st.pose.inverse().get_translation();
    return os;
}

/// The Model class implementing a constant velocity model. It also contains a damping factor
/// that will attenuate the velocities at a rate of damp/s .
/// @ingroup constantvelocitygroup
class Model {
public:
    TooN::Vector<State::STATE_DIMENSION> sigma;
    // dampening of velocity
    double damp;

    Model(void){
        TooN::Zero(sigma);
        damp = 1;
    }

    // Jacobian has pos, rot, vel, angularVel in this order
    TooN::Matrix<State::STATE_DIMENSION> getJacobian(const State & state, double dt){
            TooN::Matrix<State::STATE_DIMENSION> result;
            TooN::Identity(result);
            TooN::Identity(result.slice<0,6,6,6>(), dt);
            return result;
    }

    void updateState( State & state, const double dt ){
        // full velocity vector
        TooN::Vector<6> vel;
        vel.slice<0,3>() = state.velocity;
        vel.slice<3,3>() = state.angularVelocity;

        // update translational components
        state.pose = TooN::SE3::exp(vel * dt) * state.pose;
        // dampen velocitys
        double attenuation = pow(damp, dt);
        state.velocity *= attenuation;
        state.angularVelocity *= attenuation;
    }

    TooN::Matrix<State::STATE_DIMENSION> getNoiseCovariance( double dt ){
        TooN::Matrix<State::STATE_DIMENSION> result;
        TooN::Zero(result);
        double dt2 = dt * dt * 0.5;
        double dt3 = dt * dt * dt * 0.3333333333333;
        for(unsigned int i = 0; i < 6; i++){
            result(i,i) = dt * sigma[i] + dt3 * sigma[i+6];
            result(i+6,i) = result(i,i+6) = dt2 * sigma[i+6];
            result(i+6, i+6) = dt * sigma[i+6];
        }
        return result;
    }

    void updateFromMeasurement( State & state, const TooN::Vector<State::STATE_DIMENSION> & innovation ){
        state.pose = TooN::SE3::exp(innovation.slice<0,6>()) * state.pose;
        state.velocity = state.velocity + innovation.slice<6,3>();
        state.angularVelocity = state.angularVelocity + innovation.slice<9,3>();
    }
};

} // namespace ConstantVelocity

} // namespace tag

#endif
