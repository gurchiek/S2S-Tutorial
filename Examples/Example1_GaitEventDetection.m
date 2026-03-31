%% EXAMPLE 1 - GAIT EVENT DETECTION

% generates Figure 4 in the manuscript

% applies algorithms 2 and 3 to gait data for estimating sagittal axis of
% shank in from gait data (optionally: normal, toe in, or toe out gait).
% Uses characteristic shank sagittal angular velocity signal to sign
% correct the axes. Then detects gait events

clear
close all
clc

addpath(genpath(fullfile('..','algorithms')));
addpath(genpath(fullfile('..','helpers')));

%% INIT

% name of shank sensor to calibrate: 'IMU1' or 'IMU2'
sensor = 'IMU2';

% name of walking trial: 'Walk','WalkToeIn', or 'WalkToeOut'
walking_trial_name = 'Walk'; 

% load data
load('data.mat');

%% AXIS IDENTIFICATION FROM GAIT DATA

% gyro-based algorithm 2
gyro = estimateAxisFromUniaxialRotation(data.(walking_trial_name).imu.(sensor));

% accel-based algorithm 3
accel = estimateAxisFromPlanarMotion(data.(walking_trial_name).imu.(sensor));

%% SIGN CORRECTION

% project gyro data onto the 2 estimated axes
Yw = data.(walking_trial_name).imu.(sensor).gyro;
wz.gyro = gyro.axis.estimate' * Yw;
wz.accel = accel.axis.estimate' * Yw;

% get peaks above 1 rad/s
peaks.gyro = mean(findpeaks(wz.gyro,'MinPeakHeight',1));
peaks.accel = mean(findpeaks(wz.accel,'MinPeakHeight',1));

% get mean valley below -1 rad/s
valleys.gyro = mean(findpeaks(-wz.gyro,'MinPeakHeight',1));
valleys.accel = mean(findpeaks(-wz.accel,'MinPeakHeight',1));

% correct sign of sagittal plane gyro signal
% positive peaks should be greater than negative if sign correct
if valleys.gyro > peaks.gyro
    wz.gyro = -wz.gyro;
end
if valleys.accel > peaks.accel
    wz.accel = -wz.accel;
end

%% GAIT EVENT DETECTION

% number of strides to visualize
Nstrides = 4;

% gait event detection
t = data.(walking_trial_name).imu.time;
sf = data.(walking_trial_name).imu.samplingFrequency;
events.gyro = detectGaitEventsShankGyro(t,wz.gyro,sf);
events.accel = detectGaitEventsShankGyro(t,wz.accel,sf);

% delete first and last strides since may be incomplete
events.gyro(1) = [];
events.gyro(end) = [];
events.accel(1) = [];
events.accel(end) = [];

% foot contact times
tfc.gyro = [events.gyro.firstFootContactTimestamp];
tfc.accel = [events.accel.firstFootContactTimestamp];

% foot contact indices
ifc.gyro = zeros(1,length(tfc.gyro));
for k = 1:length(tfc.gyro); [~,ifc.gyro(k)] = min(abs(t - tfc.gyro(k))); end
ifc.accel = zeros(1,length(tfc.accel));
for k = 1:length(tfc.accel); [~,ifc.accel(k)] = min(abs(t - tfc.accel(k))); end

% foot off times
tfo.gyro = [events.gyro.firstFootOffTimestamp];
tfo.accel = [events.accel.firstFootOffTimestamp];

% foot  off indices
ifo.gyro = zeros(1,length(tfo.gyro));
for k = 1:length(tfo.gyro); [~,ifo.gyro(k)] = min(abs(t - tfo.gyro(k))); end
ifo.accel = zeros(1,length(tfo.accel));
for k = 1:length(tfo.accel); [~,ifo.accel(k)] = min(abs(t - tfo.accel(k))); end

% plot raw data
fig1 = figure;
fig1.Position = [1258 529 217 403];
axis_names = {'x','y','z'};
tstart = events.gyro(1).firstNegativeGoing0CrossingTimestamp;
tend = events.gyro(Nstrides).secondPositiveGoing0CrossingTimestamp;
duration = tend-tstart;
for k = 1:3
    sp = subplot(3,1,k);
    plot(t,Yw(k,:),'LineWidth',1,'Color',0.5*[1 1 1])
    ylabel(['Gyro-' axis_names{k} ' (rad/s)'])
    sp.Box = 'off';
    xlim([0 duration])
    if k < 3
        xticks([])
    else
        xticks(0:1:duration)
        xlabel('Time (s)')
    end

end

% gyro-based colors
transparency = 0.1;
wcolor = (1-transparency)*[0.3216, 0.1765, 0.502] + transparency*[1 1 1];

% accel-based colors
transparency = 0.4;
acolor = (1-transparency)*[0.9608, 0.4, 0] + transparency*[1 1 1];

% plot calibrated sagittal axis data
fig2 = figure;
fig2.Position =  [1258 331 217 135];
plot(t,wz.gyro,'k','LineWidth',1.5,'Color',wcolor)
hold on
plot(t,wz.accel,'LineWidth',1.5,'Color',acolor)
ylabel('Shank \omega_z (rad/s)')
xlabel('Time (s)')

% scatter plot foot contact times
scatter(tfc.gyro,wz.gyro(ifc.gyro),50,'v','filled','MarkerFaceColor',wcolor)
scatter(tfo.gyro,wz.gyro(ifo.gyro),50,'^','filled','MarkerFaceColor',wcolor)
scatter(tfc.accel,wz.accel(ifc.accel),100,'v','MarkerEdgeColor',acolor,'LineWidth',1)
scatter(tfo.accel,wz.accel(ifo.accel),100,'^','MarkerEdgeColor',acolor,'LineWidth',1)
fig2.Children.Box = 'off';
xlim([0 duration])
xticks(0:1:duration)

%% REPORT

fprintf('-Angle between axes estimated from methods 1 and 2: %4.2f degrees\n',acosd(gyro.axis.estimate'*accel.axis.estimate))
fprintf('-Method 3 (accel-based) axis estimate uncertainty relative to method 2 (gyro-based): %4.2f\n',accel.axis.uncertainty/gyro.axis.uncertainty)
fprintf('-Mean absolute error between foot contact time estimated with methods 2 and 3: %5.2f ms\n',mean(abs(tfc.accel-tfc.gyro))*1000)
fprintf('-Mean absolute error between foot off time estimated with methods 2 and 3: %5.2f ms\n',mean(abs(tfo.accel-tfo.gyro))*1000)
fprintf('-Number of strides: %d\n',length(tfc.gyro))
