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
        position = SE3();
        Zero(angularVelocity);
        Zero(velocity);
        Identity(covariance);
    }

    inline void resetVelocity(void){
        Zero(angularVelocity);
        Zero(velocity);
    }

    static const int STATE_DIMENSION = 12;
    SE3 pose;
    Vector<3> angularVelocity;
    Vector<3> velocity;
    Matrix<STATE_DIMENSION> covariance;
};

/// operator to print out instances of State in a useable manner.
/// @ingroup constantvelocitygroup
inline std::ostream & operator<< (std::ostream & os , const State & st){
    os << st.pose.ln() << st.velocity << st.angularVelocity << st.pose.inverse().get_translation();
    return os;
}

/// The Model class implementing a constant velocity model. It also contains a damping factor
/// that will attenuate the velocities at a rate of damp/s .
/// @ingroup constantvelocitygroup
class Model {
public:
    Vector<State::STATE_DIMENSION> sigma;
    // dampening of velocity
    double damp;

    Model(void){
        Zero(sigma);
        damp = 1;
    }

    // Jacobian has pos, rot, vel, angularVel in this order
    Matrix<State::STATE_DIMENSION> getJacobian(const State & state, double dt){
            Matrix<State::STATE_DIMENSION> result;
            Identity(result);
            Identity(result.slice<0,6,6,6>(), dt);
            return result;
    }

    void updateState( State & state, const double dt ){
        // full velocity vector
        Vector<6> vel;
        vel.slice<0,3>() = state.velocity;
        vel.slice<3,3>() = state.angularVelocity;

        // update translational components
        state.pose = SE3::exp(vel * dt) * state.pose;
        // dampen velocitys
        double attenuation = pow(damp, dt);
        state.velocity *= attenuation;
        state.angularVelocity *= attenuation;
    }

    Matrix<State::STATE_DIMENSION> getNoiseCovariance( double dt ){
        Matrix<State::STATE_DIMENSION> result;
        Zero(result);
        double dt2 = dt * dt * 0.5;
        double dt3 = dt * dt * dt * 0.3333333333333;
        for(unsigned int i = 0; i < 6; i++){
            result(i,i) = dt * sigma[i] + dt3 * sigma[i+6];
            result(i+6,i) = result(i,i+6) = dt2 * sigma[i+6];
            result(i+6, i+6) = dt * sigma[i+6];
        }
        return result;
    }

    void updateFromMeasurement( State & state, const Vector<State::STATE_DIMENSION> & innovation ){
        state.pose = SE3::exp(innovation.slice<0,6>()) * state.pose;
        state.velocity = state.velocity + innovation.slice<6,3>();
        state.angularVelocity = state.angularVelocity + innovation.slice<9,3>();
    }
};

} // namespace ConstantVelocity

} // namespace tag

#endif
