#include <iostream>
#include <fstream>
#include <eigen3/Eigen/Core>
#include "imu_integrator/imu_integrator.h"

#define PI 3.14159265
using namespace std;
using namespace IMUDynamics;

/** 
 * Testing IMU Integrator, the testing case is circular movement with angular speed 1 r/s (2*pi rad/s).
 * To ameliorate the accumulation of numerical error, and get high accuracy, high sample rate is used, 
 * with a sampling period of 5e-5 s. 
 * The outout trajectory is TUM format, with each line: timestamp tx ty tz qx qy qz qw .
 * The trajectory can be visualized with evo: "evo_traj tum trajectory.txt --plot" .
 */
 
int main()
{
    imuState<float> initState;
    
    initState.p_I_G.setZero();
    initState.v_I_G = Eigen::Vector3f(2*PI, 0, 0);
    initState.q_IG.setIdentity();
    
    imuReading<float> reading;
    reading.omega = Eigen::Vector3f(0, 0, 2*PI);
    reading.a = Eigen::Vector3f(0, 4*PI*PI, 0);
    reading.dT = 0.00005;
    
    IMUIntegrator integrator(initState);
    
    ofstream f("trajectory.txt");
    f<<fixed;
    
    cout<<"Testing integrator ..."<<endl;
    double current_time = 0.0;
    while(current_time < 3.0)
    {
        integrator.push_back(reading, current_time);
        
        if(integrator.isReady())
        {
            Vector3f t = integrator.getPosition();
            Vector4f q = integrator.getOrientation();
            
            f << setprecision(6) << current_time << setprecision(7) << " " << t(0) << " " << t(1) << " " << t(2)
                << " " << q[0] << " " << q[1] << " " << q[2] << " " << q[3] << endl;
        }
        
        current_time += reading.dT;
    }
    cout<<"done."<<endl;
    
    return 0;
}
