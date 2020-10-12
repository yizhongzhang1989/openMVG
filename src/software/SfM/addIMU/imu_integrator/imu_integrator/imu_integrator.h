/** 
 * This is based on the msckf_mono implementation from https://github.com/daniilidis-group/msckf_mono
 */

#pragma once
#ifndef IMU_INTEGRATOR_H_
#define IMU_INTEGRATOR_H_

#include "msckf.h"

namespace IMUDynamics
{

using namespace msckf_mono;

class IMUIntegrator
{
public:

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    /** 
     * Default constructor, the "push_back()" needs 200 readings to initialize and propagate dynamics. 
     */
    IMUIntegrator();
    
    /** 
     * Constructor with initial state, the "push_back()" can propagate dynamics immediately. 
     */
    IMUIntegrator(imuState<float>& init_imu_state_);
    
    /** 
     * Input IMU readings and propagate IMU dynamics. If the initial state is not provided, 
     * the first 200 readings are used to initialize the integrator, including estimating 
     * gravity and biases. 
     */
    void push_back(imuReading<float>& current_imu, double cur_imu_time);
    
    /** 
     * Getting IMU position. 
     */
    Eigen::Vector3f getPosition();
    
    /** 
     * Getting IMU orientation. 
     */
    Eigen::Vector4f getOrientation();
    
    /** 
     * Getting IMU full-state. 
     */
    imuState<float> getIMUState();
    
    /** 
     * Whether the initialization has finished. 
     */
    bool isReady() { return imu_calibrated_;}
    
private:
    
    /** 
     * Estimating gravity and biases, and also estimating the initial orientation. 
     */
    void initialize_imu();
    
    /** 
     * Setup the MSCKF. The IMU propagation is done by MSCKF. 
     */
    void setup_msckf(const std::string config = "");
    
    /** 
     * Flag whether the initialization has finished.
     */
    std::atomic<bool> imu_calibrated_;
    
    /** 
     * Queue of IMU readings used for IMU initialization. 
     */
    std::vector<std::tuple<double, msckf_mono::imuReading<float>>> imu_queue_;
    
    /** 
     * MSCKF, see the official document. 
     */
    MSCKF<float> msckf_;
    
    /** 
     * initial IMU state, either from the constructor providing the initial state 
     * or from online initialization (the "initialize_imu()" called from "push_back()" )
     */
    imuState<float> init_imu_state_;
    
    /** 
     * Temporary variables used for MSCKF construction (required by the constructor), not used. 
     */
    Camera<float> camera_;
    
    /** 
     * Temporary variables used for MSCKF construction (required by the constructor), not used. 
     */
    noiseParams<float> noise_params_;
    
    /** 
     * Temporary variables used for MSCKF construction (required by the constructor), not used. 
     */
    MSCKFParams<float> msckf_params_;
    
};  // class IMUIntegrator

}  // namespace msckf_mono

#endif  // IMU_INTEGRATOR_H_
