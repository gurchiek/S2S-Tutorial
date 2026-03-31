%% EXAMPLE 2: SHANK CALIBRATION USING SUPINE-SWING-STAND SEQUENCE

% generates Figure 5 in the manuscript

clear
close all
clc

addpath(genpath(fullfile('..','algorithms')));
addpath(genpath(fullfile('..','helpers')));

%% TRIAL SETTINGS

% name of shank sensor to calibrate: 'IMU1' or 'IMU2'
sensor = 'IMU2';

% walking trial in which to calculate shank rotational kinetic energy
% can be either 'Walk', 'WalkToeIn', or 'WalkToeOut'
walk_trial = 'Walk';

%% LOAD DATA AND CONSTRUCT STRUCT INPUTS TO S2S FUNCTION

% load data
load('data.mat');

% construct supine struct for shank s2s function
supine.accel = data.Supine.imu.(sensor).accel;
supine.gyro = data.Supine.imu.(sensor).gyro;
supine.sf = data.Supine.imu.samplingFrequency;

% construct swing struct for shank s2s function
swing.accel = data.SeatedLegSwing.imu.(sensor).accel;
swing.gyro = data.SeatedLegSwing.imu.(sensor).gyro;
swing.sf = data.SeatedLegSwing.imu.samplingFrequency;

% construct su struct for shank s2s function
stand.accel = data.Standing.imu.(sensor).accel;
stand.gyro = data.Standing.imu.(sensor).gyro;
stand.sf = data.Standing.imu.samplingFrequency;

%% CALIBRATE

results = s2sShankCalibrationSupineSwingStand(supine,swing,stand);

%% GET TRUE SENSOR TO SEGMENT ORIENTATION

% use data during HipStar trial for ground truth alignment
timu = data.HipStar.imu.time;
wimu = data.HipStar.imu.(sensor).gyro - results.gyro.bias;

% reconstruct shank motion from marker data
omc = reconstructShankMotionFromMarkerData(data.Standing.trc,data.HipStar.trc);
womc = omc.wB;
tomc = omc.time;

% low pass filter at 4 hz
wimu = bwfilt(wimu,4,data.HipStar.imu.samplingFrequency,'low',4);
womc = bwfilt(womc,4,omc.sf,'low',4);

% time-synchronize omc and imu based on angular velocity magnitude
sync = digitalSignalSynchronizationRMS(vecnorm(wimu),timu,vecnorm(womc),tomc,0.001,1,-inf,inf);
tomc = sync.t2;

% interpolate data to same time array
t = max([tomc(1),timu(1)]):1/225:min([tomc(end),timu(end)]);
womc = interp1(tomc,womc',t,'pchip')';
wimu = interp1(timu,wimu',t,'pchip')';

% use SVD soln to Wahba problem to estimate true S2S orientation with
% product of magnitudes as weight
Usensor = normalize(wimu,1,'norm');
Ubody = normalize(womc,1,'norm');
w = vecnorm(wimu) .* vecnorm(womc);

% solve wahbas problem
omcWahba = estimateOrientationFromAxisEstimates(Ubody,Usensor,w);

%% ERROR STATS

% get omc and imu-based s2s orientations
Romc = omcWahba.R;
Rimu = results.R_S2B;

% error angle between Rest and Rtrue
Rerror = Rimu' * Romc;
error_angle = acosd((trace(Rerror) - 1)/2);
fprintf('\n-Angle of error rotation matrix (OMC vs. IMU) = %5.3f degrees\n',error_angle)

% angle between true and estimated x-, y-, and z-axes
fprintf('-Error angle between estimated and true x-axis: %4.2f degrees\n',acosd(Romc(1,:)*Rimu(1,:)'))
fprintf('-Error angle between estimated and true y-axis: %4.2f degrees\n',acosd(Romc(2,:)*Rimu(2,:)'))
fprintf('-Error angle between estimated and true z-axis: %4.2f degrees\n',acosd(Romc(3,:)*Rimu(3,:)'))

%% SHANK ROTATIONAL KINETIC ENERGY DURING GAIT

% this illustrates the effect of sensor to segment orientation for a
% biomechanical outcome that combines all three sense axes of the
% gyroscope: rotational kinetic energy

% get gyro data during gait
walk.imu.time = data.(walk_trial).imu.time;
walk.imu.gyro = data.(walk_trial).imu.(sensor).gyro;
walk.imu.sf = data.(walk_trial).imu.samplingFrequency;

% calibrate gyro
walk.imu.wB = Rimu * (walk.imu.gyro - results.gyro.bias);

% reconstruct shank motion from marker data
walk.omc = reconstructShankMotionFromMarkerData(data.Standing.trc,data.(walk_trial).trc);

% low pass filter
walk.imu.gyro = bwfilt(walk.imu.gyro,6,walk.imu.sf,'low',4);
walk.imu.wB = bwfilt(walk.imu.wB,6,walk.imu.sf,'low',4);
walk.omc.w = bwfilt(walk.omc.wB,6,walk.omc.sf,'low',4);

% inertia matrix for shank (taken from OpenSim Rajagopal model)
I = diag([0.0504, 0.0051, 0.0511]);

% calculate rotational kinetic energy based on calibrated and uncalibrated
% gyro data and omc-based angular velocity
Tuncal = dot(walk.imu.gyro,I * walk.imu.gyro) / 2;
Tcal = dot(walk.imu.wB,I * walk.imu.wB) / 2;
Tomc = dot(walk.omc.wB,I * walk.omc.wB) / 2;

% synchronize Tomc with Tcal
sync = digitalSignalSynchronizationRMS(Tcal/max(Tcal),walk.imu.time,Tomc/max(Tomc),walk.omc.time,0.001,1,-2,2);
walk.omc.time = sync.t2;

%% VISUALIZE AND COMPARE UNCALIBRATED VS CALIBRATED GYRO DATA

% number of strides to visualize
Nstrides = 4;

% line color and width
color.uncal = 0.7*[1 1 1];
color.cal = [0 0 0];
color.omc = 0.9*[0.2745, 0.5059, 0.5059] + (1-0.9)*[1 1 1];
lineWidth = 2;
lineStyle.cal = '-';
lineStyle.uncal = '-';
lineStyle.omc = '-';

% get gait events
events = detectGaitEventsShankGyro(walk.imu.time,walk.imu.wB,walk.imu.sf);
events(1) = []; % delete first stride since might be incomplete

% get time limits consistent with requested num strides to plot
tstart = events(1).firstNegativeGoing0CrossingTimestamp;
tend = events(Nstrides).secondPositiveGoing0CrossingTimestamp;
duration = tend-tstart;

% visualize x, y, and z axes ang vel from omc, and cal/uncal imu
fig1 = figure;
fig1.Position = [1366 507 217 537];
axis_names = {'x','y','z'};
for k = 1:3

    % uncal and omc
    sp = subplot(4,1,k);
    plot(walk.imu.time-tstart,walk.imu.wB(k,:),'LineWidth',lineWidth,'Color',color.cal,'LineStyle',lineStyle.cal)
    hold on
    plot(walk.imu.time-tstart,walk.imu.gyro(k,:),'LineWidth',lineWidth,'Color',color.uncal,'LineStyle',lineStyle.uncal)
    plot(walk.omc.time-tstart,walk.omc.wB(k,:),'LineWidth',lineWidth,'Color',color.omc,'LineStyle',lineStyle.omc)
    ylabel(['\omega ' axis_names{k} '-axis (rad/s)'])
    xlim([0 duration])
    xticks([])
    sp.Box = 'off';

end

% rotational KE
sp = subplot(4,1,4);
plot(walk.imu.time-tstart,Tcal,'LineWidth',lineWidth,'Color',color.cal,'LineStyle',lineStyle.cal)
hold on
plot(walk.imu.time-tstart,Tuncal,'LineWidth',lineWidth,'Color',color.uncal,'LineStyle',lineStyle.uncal)
plot(walk.omc.time-tstart,Tomc,'LineWidth',lineWidth,'Color',color.omc,'LineStyle',lineStyle.omc)
xlabel('Time (s)')
ylabel('Rot. KE (J)')
xlim([0 duration])
xticks(0:1:duration)
sp.Box = 'off';
