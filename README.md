# S2S-Tutorial
Example code that accompanies the tutorial on estimating sensor-to-segment orientation published in the *Journal of Biomechanics* - please cite the following: TBD

Dependencies: Matlab + Signal Processing Toolbox

# Examples
1. Gait Event Detection - Implements data-driven example 1 in the paper and generates Figure 4
2. Shank Calibration - Implements data-driven example 2 in the paper and generates Figure 5
3. Knee Flexion Axis Estimation - Implements application of algorithms 2 and 3 to example calibration motion 1 in the paper and generates Figures 2 and 3
4. Shank Calibration Simulation - Implements the shank sensor-to-segment orientation estimation portion of data-driven example 2 in the paper applied to simulated data with user-specified sensor characteristics (noise, bias, etc.)

# Algorithms
1. estimateAxisFromStaticPosture - Implements algorithm 1 in the paper
2. estimateAxisFromUniaxialRotation - Implements algorithm 2 in the paper
3. estimateAxisFromPlanarMotion - Implements algorithm 3 in the paper
4. s2sShankCalibrationSupineSwingStand - Implements algorithm 4 in the paper

# Data
IMU and Marker data are stored in the Matlab data file 'data.mat'. It is a struct organized as per:

'''text
data.(trial_name).imu.time = 1xN array of timestamps
                     .numSamples
                     .samplingFrequency
                     .(imu_name).gyro = 3xN array of gyroscope data (rad/s)
                                .accel = 3xN array of accelerometer data (m/s/s)
                  .trc.numMarkers
                      .numSamples
                      .time = 1xN array of timestamps
                      .samplingFrequency
                      .marker.(marker_name).position = 3xN array of marker position data (m)

trial_names: Standing, Supine, SeatedLegSwing, Walk, WalkToeIn, WalkToeOut, HipStar (see paper for details)

imu_names: IMU1, IMU2, IMU3 (see paper for details)

marker_names: IMU1, IMU2, IMU3, LMAL, FIBHEAD, MMAL, TIBTUB, SHANK1, SHANK2, SHANK3, LEPI, MEPI, THIGH1, THIGH2, THIGH3, THIGH4 (see paper for details)
                               
