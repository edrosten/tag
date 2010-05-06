#ifndef TAG_MEASUREMENTS_H
#define TAG_MEASUREMENTS_H

#include <TooN/se3.h>

namespace tag {

/// @defgroup measurementsgroup Measurement classes
/// This group contains model independend measurement classes that can be used with a
/// variety of models such as ConstantPosition or ConstantVelocity. They typically
/// need to be parameterized on the type of the state.
/// @ingroup kalmanfiltergroup

/// An incremental measurement consisting of the 6-vector parameterizing a correction to be
/// left multiplied to the current pose. This is typically the result of an optimization run.
/// @ingroup measurementsgroup
template <class State>
class IncrementalPose {
public:
    static const int M_DIMENSION = 6;

    TooN::Matrix<M_DIMENSION> covariance;
    TooN::Matrix<M_DIMENSION,State::STATE_DIMENSION> jacobian;
    TooN::Vector<M_DIMENSION> measurement;

    IncrementalPose(void){
        covariance = TooN::Identity;
        jacobian = TooN::Zeros;
        jacobian.template slice<0,0,6,6>() = TooN::Identity;
        measurement = TooN::Zeros;
    }

    TooN::Matrix<M_DIMENSION,State::STATE_DIMENSION> & getMeasurementJacobian( const State & state ){
        return jacobian;
    }

    TooN::Matrix<M_DIMENSION> & getMeasurementCovariance( const State & state ){
        return covariance;
    }

    TooN::Vector<M_DIMENSION> & getInnovation( const State & state ){
        return measurement;
    }

    void setMeasurement( const TooN::SE3<> & increment ){
        measurement = increment.ln();
    }
};

/// An absolute pose measurement in world coordinates. This is transformed
/// into an incremental motion given the state of the filter.
/// @ingroup measurementsgroup
template <class State>
class WorldPose {
public:
    static const int M_DIMENSION = 6;

    TooN::Matrix<M_DIMENSION> covariance;
    TooN::Matrix<M_DIMENSION,State::STATE_DIMENSION> jacobian;
    TooN::SE3<> measurement;

    WorldPose(void){
        covariance = TooN::Identity;
        jacobian = TooN::Zeros;
        jacobian.template slice<0,0,6,6>() = TooN::Identity;
    }

    TooN::Matrix<M_DIMENSION,State::STATE_DIMENSION> & getMeasurementJacobian( const State & state ){
        return jacobian;
    }

    const TooN::Matrix<M_DIMENSION> getMeasurementCovariance( const State & state ) const {
        TooN::Matrix<M_DIMENSION> localCovariance = covariance;
        state.pose.adjoint(localCovariance);
        return localCovariance;
    }

    const TooN::Vector<M_DIMENSION> getInnovation( const State & state ) const {
        return (measurement * state.pose.inverse()).ln();
    }
};

/// An absolute position (= translation) measurement in world coordinates. This is transformed
/// into an incremental motion given the state of the filter. Possible sensors include GPS.
/// @ingroup measurementsgroup
template <class State>
class WorldPosition {
public:

    static const int M_DIMENSION = 3;
    TooN::Vector<M_DIMENSION> position;
    TooN::Matrix<M_DIMENSION> covariance;
    TooN::Matrix<M_DIMENSION,State::STATE_DIMENSION> jacobian;

    WorldPosition(void){
        covariance = TooN::Identity;
        position = TooN::Zeros;
        jacobian = TooN::Zeros;
        jacobian.template slice<0,0,3,3>() = TooN::Identity;
    }

    const TooN::Matrix<M_DIMENSION,State::STATE_DIMENSION> & getMeasurementJacobian( const State & state ) const {
        return jacobian;
    }

    const TooN::Matrix<M_DIMENSION> getMeasurementCovariance( const State & state ) const {
        TooN::Matrix<2 * M_DIMENSION> localCovariance = TooN::Zeros;
        localCovariance.template slice<0,0,3,3>() = covariance;
        state.pose.adjoint(localCovariance);
        return localCovariance.template slice<0,0,3,3>();
    }

    const TooN::Vector<M_DIMENSION> getInnovation( const State & state ) const {
        /// the negative vector corresponds to the left transformation to get from the current reference
        /// frame to the new reference frame where position is the origin.
        return -(state.pose * position);
    }

    void setCovariance( double sigma ){
        covariance = TooN::Identity * sigma;
    }

    void setCovariance( const TooN::Vector<M_DIMENSION> & sigma ){
        covariance = TooN::Zeros;
        for(int i = 0; i < M_DIMENSION; ++i){
            covariance(i,i) = sigma[i];
        }
    }
};

/// A angular velocity measurement in local coordinates. Possible sensors include gyroscopes.
/// @ingroup measurementsgroup
template <class State>
class AngularVelocity {
public:

    static const int M_DIMENSION = 3;

    TooN::Vector<3> gyro;
    TooN::Matrix<M_DIMENSION> covariance;
    TooN::Matrix<M_DIMENSION,State::STATE_DIMENSION> jacobian;

    AngularVelocity(void){
        covariance = TooN::Identity;
        gyro = TooN::Zeros;
        // angularVelocity jacobian
        jacobian = TooN::Zeros;
        jacobian.template slice<0,9,3,3>() = TooN::Identity;
    }

    TooN::Matrix<M_DIMENSION,State::STATE_DIMENSION> & getMeasurementJacobian( const State & state ){
        return jacobian;
    }

    TooN::Matrix<M_DIMENSION> & getMeasurementCovariance( const State & state ){
        return covariance;
    }

    TooN::Vector<M_DIMENSION> getInnovation( const State & state ){
        // angular velocity
        /// @todo think about the -gyro and either document or remove
        return  (-gyro - state.angularVelocity);
    }

    void setCovariance( double sigma ){
        covariance = TooN::Identity * sigma;
    }
};

/// Measuring rotation by comparing a measured direction to a known direction in world space.
/// Possible sensors include gravity or magnetic field direction.
/// @ingroup measurementsgroup
template <class State>
class WorldDirection {
public:

    static const int M_DIMENSION = 3;

    TooN::Vector<3> measurement;
    TooN::Vector<3> reference;

    TooN::Matrix<M_DIMENSION> covariance;

    WorldDirection(void){
        covariance = TooN::Identity;
        measurement = TooN::Zeros;
        reference = TooN::Zeros;
    }

    TooN::Matrix<M_DIMENSION,State::STATE_DIMENSION> getMeasurementJacobian( const State & state ){
        TooN::Matrix<M_DIMENSION,State::STATE_DIMENSION> result = TooN::Zeros;
        // direction jacobian
        TooN::Vector<M_DIMENSION> local = state.pose.get_rotation() * reference;
        result.template slice<0,3,3,1>() = TooN::SO3<>::generator_field(0, local ).as_col();
        result.template slice<0,4,3,1>() = TooN::SO3<>::generator_field(1, local ).as_col();
        result.template slice<0,5,3,1>() = TooN::SO3<>::generator_field(2, local ).as_col();
        return result;
    }

    TooN::Matrix<M_DIMENSION> & getMeasurementCovariance( const State & state ){
        return covariance;
    }

    TooN::Vector<M_DIMENSION> getInnovation( const State & state ){
        return (measurement - (state.pose.get_rotation() * reference));
    }

    void setCovariance( double sigma ){
        covariance = TooN::Identity * sigma;
    }
};

} // namespace tag

#endif
