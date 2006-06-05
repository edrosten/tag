#ifndef TAG_FILTERKALMANFILTER_H
#define TAG_FILTERKALMANFILTER_H

#include <TooN/TooN.h>
#include <TooN/helpers.h>
#include <TooN/LU.h>

namespace tag {

/**
@defgroup kalmanfiltergroup Kalman Filter
A basic Kalman filter implementation and various state, process models and measurement functions.

Template class providing a basic implementation of the Kalman filter.
The state and the process model are both template parameter classes
to keep it flexible. Both parameters have to implement a certain interface to make the filter work.
@code
class State {
public:
    const static int STATE_DIMENSION =  ??;  // dimension of filter state
    TooN::Matrix<STATE_DIMENSION> covariance;      // covariance of filter state
};

class Model {
public:
    // return process model jacobian for a given state and time delta (typically A)
    TooN::Matrix<State::STATE_DIMENSION> getJacobian(const State & state, double dt);
    // update the state for a given time delta (not all states are actually x = Ax, therefore this is a function)
    void updateState( State & state, const double dt );
    // return process noise matrix for given time delta (typically Q)
    TooN::Matrix<State::STATE_DIMENSION> getNoiseCovariance( double dt );
    // update the state from an innovation. the innovation was computed by the filter based on measurement etc.
    void updateFromMeasurement( State & state, const TooN::Vector<State::STATE_DIMENSION> & innovation );
};
@endcode
Measurements are incorporated through the template member function
template<class Measurement> void KalmanFilter<class State, class Model>::filter(Measurement & m);
where class Measurement also has to implement a certain protocol:
@code
class Measurement {
public:
    static const int M_DIMENSION =  ??;  // dimension of measurement
    // return measurement jacobian, from state -> measurement
    Matrix<M_DIMENSION,State::STATE_DIMENSION> getMeasurementJacobian( const State & state );
    // return measurement noise covariance
    Matrix<M_DIMENSION> getMeasurementCovariance( const State & state );
    // return the innovation, the difference between actual measurement and the measurement prediction based on the state
    Vector<M_DIMENSION> getInnovation( const State & state );
};
@endcode
All of the member functions take the state as parameters, because the returned values are typically
functions of the state in some form.

Basically, the three classes State, Model, Measurement have to know about each other and work together.
However, splitting them apart allows one to change models and use multiple measurement functions for a
single kind of state. That simplifies sensor fusion and SCAAT style use of the Kalman filter.

The following example demonstrates how to use the filter classes.
@code
KalmanFilter<ConstantVelocity::State, ConstantVelocity::Model> filter;
filter.state.pose = // Initial pose
filter.state.covariance = // Initial covariance

while(true){
    double deltaT = 0.1; // interval between measurements
    filter.predict( deltaT );

    ConstantVelocity::CameraMeasurement m;
    m.measurement = // update vector = ln() of the correction SE3
    m.covariance =  // your measurement covariance
    filter.filter(m);
}
@endcode
*/

/**
the basic template class implementing the Kalman Filter, see @ref kalmanfiltergroup documentation for details.
@ingroup kalmanfiltergroup
*/
template<class State, class Model>
class KalmanFilter{
public:

    KalmanFilter(){
        TooN::Identity(identity);
    }

    /// predicts the state by applying the process model over the time interval dt
    /// @param[in] dt time interval
    void predict(double dt){
        TooN::Matrix<State::STATE_DIMENSION> A = model.getJacobian( state, dt );
        model.updateState( state, dt );
        state.covariance = TooN::transformCovariance(A, state.covariance) + model.getNoiseCovariance( dt );
        TooN::Symmetrize(state.covariance);
    }

    /// incorporates a measurement
    /// @param[in] m the measurement to add to the filter state
    template<class Measurement> void filter(Measurement & m){
        TooN::Matrix<Measurement::M_DIMENSION,State::STATE_DIMENSION> H = m.getMeasurementJacobian( state );
        TooN::Matrix<Measurement::M_DIMENSION> R = m.getMeasurementCovariance( state );
        TooN::Matrix<Measurement::M_DIMENSION> I = TooN::transformCovariance(H, state.covariance) + R;
        TooN::LU<Measurement::M_DIMENSION> lu(I);
        TooN::Matrix<State::STATE_DIMENSION, Measurement::M_DIMENSION> K = state.covariance * H.T() * lu.get_inverse();
        TooN::Vector<Measurement::M_DIMENSION> innovation = m.getInnovation( state );
        TooN::Vector<State::STATE_DIMENSION> stateInnovation = K * innovation;
        model.updateFromMeasurement( state, stateInnovation );
        state.covariance = (identity - K * H) * state.covariance;
        Symmetrize( state.covariance );
    }

    /// identity matrix of the right size, used in the measurement equations
    TooN::Matrix<State::STATE_DIMENSION> identity;
    /// the current state of the filter
    State state;
    /// the process model used by the filter
    Model model;
};

} // namespace tag

#endif
