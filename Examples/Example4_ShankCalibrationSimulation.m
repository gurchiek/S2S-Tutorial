clear
close all
clc

addpath(genpath(fullfile('..','algorithms')));
addpath(genpath(fullfile('..','helpers')));

%% SETTIINGS

% sampling rate
sf = 225; 

% clemson elevation and latitude => local gravity
elev = 221; % meters
lat = 34.68; % degrees
if license('test','Aerospace_Toolbox')
    g = gravitywgs84(elev,lat);
else
    g = 9.8066499999999994; % opensim
end

% set noise density
% Axivity AX6 uses Bosch BM160
% accel noise density = 180 - 300 micro-g's / sqrt(Hz)
% accel noise sigma = 1.8 milli-g's RMS
% gyro noise density = 0.007 deg/s / sqrt(Hz)
% gyro noise sigma = 0.07 deg/s RMS
accel_noise_density = 200e-6 * g;
gyro_noise_density = 0.007 * pi/180;
accel_noise_sigma = 1.8e-3 * g;
gyro_noise_sigma = 0.07 * pi/180;

% set noise standard deviation
% either asig = accel_noise_sigma or = accel_noise_density * sqrt(sf/2)
% approximate noise variance from static trials:
%   accel ~ 0.003 (m/s/s)^2
%   gyro ~ 0.3e-4 (rad/s)^2
% note the reported approximate noise variance from static trials are most
% similar to accel/gyro_noise_simgma's above
asig = sqrt(0.003);
wsig = sqrt(0.3e-4);

% set accel bias
% Bosch BM160 zero-g offset is +/- 40 milli-g once on-board
baccel_range = 50 * 1e-3 * g;
baccel = (2*rand(3,1) - 1) * baccel_range;

% set gyro bias
% Bosch BM160 zero-rate offset is +/- 3 deg/s once on-board
bgyro_range = 3 * pi/180;
bgyro = (2*rand(3,1) - 1) * bgyro_range;

% set orientation of shank z axis in world frame during the swing trial
% if simulating calibration motion 1, then set to [0 0 1]'
zWswing = [1 0 1]';

% set location of sensor (A) on shank in shank's body frame taking knee
% joint center as origin. Example: rOA_B = [ 0 -0.15 0.05] would model the
% sensor is located 15 cm below and 5 cm to the right of the knee JC
rOA_B = [0 -0.15 0.05]';

% set segment-to-sensor rotation matrix s.t. vS = R * vB
% the transpose of this matrix is what we're trying to estimate
% the columns of RBS are the axes estimated using methods 1-3
Raxis = normalize([1;1;1],1,'norm');
Rangle = 30 * pi/180;
RBS = eye(3) + sin(Rangle) * skew(Raxis) + (1 - cos(Rangle)) * skew(Raxis) * skew(Raxis);

%% REPORT SETTINGS

fprintf('Simulation Settings:\n')
fprintf('\n     -Hardware:\n')
fprintf('           -Gyro noise variance = %.3e rad/s\n',wsig^2);
fprintf('           -Gyro bias = %5.4f, %5.4f, %5.4f rad/s\n',bgyro)
fprintf('                -Magnitude = %5.4f rad/s\n',vecnorm(bgyro))
fprintf('           -Accel noise variance = %.3e m/s/s\n',asig^2)
fprintf('           -Accel bias = %5.4f, %5.4f, %5.4f m/s/s\n',baccel)
fprintf('                -Magnitude = %5.4f m/s/s\n',vecnorm(baccel))
fprintf('\n     -Sensor-to-segment:\n')
fprintf('           -Representation of shank x-axis in sensor frame: %5.4f, %5.4f, %5.4f\n',RBS(:,1))
fprintf('           -Representation of shank y-axis in sensor frame: %5.4f, %5.4f, %5.4f\n',RBS(:,2))
fprintf('           -Representation of shank z-axis in sensor frame: %5.4f, %5.4f, %5.4f\n',RBS(:,3))
fprintf('           -Location of sensor relative to knee joint center in shank frame: %4.2f, %4.2f, %4.2f cm\n',rOA_B*100)
fprintf('\n     -Calibration motion settings:\n')
fprintf('           -Representation of axis of rotation in world frame during leg swing task (i.e., shank z-axis): %5.4f, %5.4f, %5.4f\n',normalize(zWswing,1,'norm'))

%% SIMULATE STATIC POSTURE 1: STANDING

% num samples
N1 = 2*sf;

% simulate sensor noise
anoise1 = asig * randn(3,N1);
wnoise1 = wsig * randn(3,N1);

% shank y axis in sensor frame
yBS = RBS(:,2);

% simulate accelerometer gyro signal in standing posture
stand.accel = g * yBS + baccel + anoise1;
stand.gyro = bgyro + wnoise1;

%% SIMULATE STATIC POSTURE 2: SUPINE

% num samples
N2 = 2*sf;

% simulate accelerometer noise
anoise2 = asig * randn(3,N2);
wnoise2 = wsig * randn(3,N1);

% shank x axis in sensor frame
xBS = RBS(:,1);

% simulate accelerometer and gyro signal in standing posture
supine.accel = g * xBS + baccel + anoise2;
supine.gyro = bgyro + wnoise2;

%% SIMULATE CALIBRATION MOTION 1

% simulate shank hanging off end of table => leg swinging about zWsing (set
% above) up to 45 deg back to -45 deg and back to 0 degrees according to a
% sinusoid with period T = 1 sec for 5 swings (5 seconds)
t = 0:1/sf:5;
N3 = length(t);
T = 1/2/pi;
amp = 45*pi/180;
angle = amp * sin(t/T);

% get rotation matrix Rswing0 s.t. vW = Rswing0 * vB 
% should have zWswing = Rswing0 * zB
zB = [0; 0; 1];
I3 = eye(3);
zWswing = normalize(zWswing,1,'norm');
temp = cross(zB,zWswing);
phi = asin(vecnorm(temp));
if phi < eps
    Rswing0 = I3;
else
    n = temp / vecnorm(temp);
    Rswing0 = I3 + sin(phi) * skew(n) + (1 - cos(phi)) * skew(n) * skew(n);
end

% get rotation matrix RB that maps vW = RBW * vB
RBW = zeros(3,3,N3);
for k = 1:N3
    RBW(:,:,k) = Rswing0 * (I3 + sin(angle(k)) * skew(zB) + (1 - cos(angle(k))) * skew(zB) * skew(zB));
end

% should have RBW(:,:,k) * [0 0 1]' = zWswing set above

% get angular rotation speeds
s = amp/T * cos(t/T);

% shank z-axis in sensor frame
zBS = RBS(:,3);

% get angular velocity in sensor frame
w3 = zeros(3,N3);
for k = 1:N3
    w3(:,k) = s(k) * zBS;
end

% get angular rotation acceleration
sdot = -amp/T/T * sin(t/T);

% get angular acceleration in sensor frame
w3dot = zeros(3,N3);
for k = 1:N3
    w3dot(:,k) = sdot(k) * zBS;
end

% get matrix KBS
KBS = zeros(3,3,N3);
for k = 1:N3
    KBS(:,:,k) = skew(w3dot(:,k)) + skew(w3(:,k)) * skew(w3(:,k));
end

% get location of sensor in sensor frame
rOA_S = RBS * rOA_B;

% simulate IMU noise
wnoise3 = wsig * randn(3,N3);
anoise3 = asig * randn(3,N3);

% simulate specific force
a3 = zeros(3,N3);
for k = 1:N3
    a3(:,k) = KBS(:,:,k) * rOA_S + g * RBS * RBW(:,:,k)' * [0; 1; 0];
end

% simulate accelerometer and gyroscope signal
swing.accel = a3 + baccel + anoise3;
swing.gyro = w3 + bgyro + wnoise3;

%% CALIBRATE

% add sampling frequency to inputs structs to calibration function
stand.sf = sf;
swing.sf = sf;
supine.sf = sf;

% calibrate
results = s2sShankCalibrationSupineSwingStand(supine,swing,stand);

%% COMPARE ESTIMATED AXES TO SIMULATED GROUND TRUTH

fprintf('\nComparison of estimated axes to simulated ground truth:\n')

fprintf('\n-Axis estimation error angle from method 1 (compared to simulation ground truth):\n')
xu_stand = results.stand.calibration.axis.estimate;
phi_stand = acosd(yBS' * xu_stand);
fprintf('     -Stand: %4.2f deg\n',phi_stand)
xu_supine = results.supine.calibration.axis.estimate;
phi_supine = acosd(xBS' * xu_supine);
fprintf('     -Supine: %4.2f deg\n',phi_supine)

% estimate z-axis using method 3
method3 = estimateAxisFromPlanarMotion(swing,struct('initial_guess',results.R_S2B(3,:)'));

fprintf('\n-Comparison of methods 2 and 3 for estimating z-axis (compared to simulation ground truth):\n')
fprintf('     -Method 2: error angle = %4.2f deg\n',acosd(results.swing.calibration.axis.estimate' * RBS(:,3)))
fprintf('           -Axis estimate: %5.4f, %5.4f, %5.4f\n',results.swing.calibration.axis.estimate)
fprintf('           -Uncertainty: %.3e\n',results.swing.calibration.axis.uncertainty)
fprintf('     -Method 3: error angle = %4.2f deg\n',acosd(method3.axis.estimate' * RBS(:,3)))
fprintf('           -Axis estimate: %5.4f, %5.4f, %5.4f\n',method3.axis.estimate)
fprintf('           -Uncertainty: %.3e\n',method3.axis.uncertainty)

% estimate z-axis using method 3 with abnormally large accel bias to
% demonstrate independence of the method on accel bias
abnormally_large_baccel = 100;
swingBias.accel = swing.accel + abnormally_large_baccel * ones(3,N3);
method3Bias = estimateAxisFromPlanarMotion(swingBias,struct('initial_guess',results.R_S2B(3,:)'));

fprintf('\n-Comparison of method 3 with and without large bias for estimating z-axis (should see no difference): bias = %6.3f m/s/s (per axis)\n',abnormally_large_baccel)
fprintf('     -Method 3 small bias: %5.4f, %5.4f, %5.4f\n',method3.axis.estimate)
fprintf('           -Uncertainty: %.3e\n',method3.axis.uncertainty)
fprintf('     -Method 3 large bias: %5.4f, %5.4f, %5.4f\n',method3Bias.axis.estimate)
fprintf('           -Uncertainty: %.3e\n',method3Bias.axis.uncertainty)
fprintf('\n---------------------------------------------------\n')



%% COMPARE ESTIMATED ORIENTATION TO SIMULATION GROUND TRUTH

fprintf('\nComparison of estimated orientation to simulated ground truth:\n')

% get true and imu-based s2s orientations
Rtrue = RBS';
Rimu = results.R_S2B;

% error angle between Rest and Rtrue
Rerror = Rimu' * Rtrue;
error_angle = acosd((trace(Rerror) - 1)/2);
fprintf('\n-Angle of error rotation matrix (simulation ground truth vs. estimated) = %4.2f degrees\n',error_angle)

% angle between true and estimated x-, y-, and z-axes
fprintf('-Error angle between estimated and simulated true x-axis: %4.2f degrees\n',acosd(Rtrue(1,:)*Rimu(1,:)'))
fprintf('-Error angle between estimated and simulated true y-axis: %4.2f degrees\n',acosd(Rtrue(2,:)*Rimu(2,:)'))
fprintf('-Error angle between estimated and simulated true z-axis: %4.2f degrees\n',acosd(Rtrue(3,:)*Rimu(3,:)'))
