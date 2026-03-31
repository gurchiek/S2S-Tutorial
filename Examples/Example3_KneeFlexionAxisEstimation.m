%% Example 3: Estimate Knee Flexion Axis

% generates Figures 2 and 3 in the manuscript

% applies algorithm 2 to gyroscope data and algorithm 3 to accelerometer
% data during a leg swing task (example calibration motion 1)

% axis estimates are compared and visualized

clear
close all
clc

addpath(genpath(fullfile('..','algorithms')));
addpath(genpath(fullfile('..','helpers')));

%% SETTINGS

% will plot only every nth point
nth = 50;

% properties of scatter plot and arrow illustrating axis
arrowWidth = 1;
arrowColor = [1 0 0];
arrowScale = 0.9;
mkrEdgeColor = [0 0 0];
mkrFaceColor = 'none';
mkrFaceAlpha = 0.03;
mkrSize = 50;
mkrLineWidth = 2;
mkrEdgeAlpha = 0.1;
limScale = 1.25;

%% INIT

% load data
load('data.mat');

% name of shank sensor to calibrate: 'IMU1' or 'IMU2'
sensor = 'IMU2';

% construct swing struct for shank s2s function
swing.accel = data.SeatedLegSwing.imu.(sensor).accel;
swing.gyro = data.SeatedLegSwing.imu.(sensor).gyro;
swing.sf = data.SeatedLegSwing.imu.samplingFrequency;


%% GET TRUE SENSOR TO SEGMENT ORIENTATION

% use data during HipStar trial for ground truth alignment
timu = data.HipStar.imu.time;
wimu = data.HipStar.imu.(sensor).gyro;

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
weight = vecnorm(wimu) .* vecnorm(womc);

% solve wahbas problem
omcWahba = estimateOrientationFromAxisEstimates(Ubody,Usensor,weight);
R = omcWahba.R;

%% LEG SWING CALIBRATION METHODS

% initial guess of shank z axis in sensor frame (for sign correction only)
options.initial_guess = R(3,:)';

% method 2 - gyro based
method2 = estimateAxisFromUniaxialRotation(swing,options);

% method 3 - accel based
method3 = estimateAxisFromPlanarMotion(swing,options);

%% COMPARE

% R from omc above maps sensor to body => its 3rd row is shank z in sensor
fprintf('-Method 2 estimate (gyro based): %4.3f, %4.3f, %4.3f\n',method2.axis.estimate)
fprintf('     -Estimate uncertainty = %.3e\n',method2.axis.uncertainty)
fprintf('-Method 3 estimate (accel based): %4.3f, %4.3f, %4.3f\n',method3.axis.estimate)
fprintf('     -Estimate uncertainty = %.3e\n',method3.axis.uncertainty)
fprintf('-Angle between axis estimates = %4.2f deg\n',acosd(method2.axis.estimate'*method3.axis.estimate))

%% VISUALIZE AXIS ESTIMATION FOR METHOD 2

% get gyro data
Y2 = swing.gyro;

% plot only every nth point
i2 = 1:5:size(Y2,2);
Y2 = Y2(:,i2);

% get estimated z axis from method 2
z2 = method2.axis.estimate;

% rotate from sensor to body
Y2 = R*Y2;
z2 = R*z2;

% set limits for visualization
mx = max(max(abs(Y2)));
lims = limScale*mx*[-1 1];

% stretch z for visualization
z2 = arrowScale*mx*z2;

% init figure
fig1 = figure;
fig1.Position = [1243 83 324 859];

% image scatter plot from positive x-axis: z-axis to left, y-axis up
sp1 = subplot(3,1,1);
scatter3(Y2(1,:),Y2(2,:),Y2(3,:),mkrSize,'filled','MarkerEdgeColor',mkrEdgeColor,'MarkerFaceColor',mkrFaceColor,'MarkerFaceAlpha',mkrFaceAlpha,...
    'LineWidth',mkrLineWidth,'MarkerEdgeAlpha',mkrEdgeAlpha)
grid on
hold on
quiver3(1,0,0,1+z2(1),z2(2),z2(3),'Color',arrowColor,'LineWidth',arrowWidth,'MaxHeadSize',0.3);
xlim(lims)
ylim(lims)
zlim(lims)
view([1 0 0])
camup([0 1 0])
zticks(0)
xticks(0)
yticks(0)
xticklabels([])
zticklabels([])
yticklabels([])
sp1.Box = 'on';

% image scatter plot from positive y-axis: x-axis up
sp2 = subplot(3,1,2);
scatter3(Y2(1,:),Y2(2,:),Y2(3,:),mkrSize,'filled','MarkerEdgeColor',mkrEdgeColor,'MarkerFaceColor',mkrFaceColor,'MarkerFaceAlpha',mkrFaceAlpha,...
    'LineWidth',mkrLineWidth,'MarkerEdgeAlpha',mkrEdgeAlpha)
grid on
hold on
quiver3(0,1,0,z2(1),1+z2(2),z2(3),'Color',arrowColor,'LineWidth',arrowWidth,'MaxHeadSize',0.3);
xlim(lims)
ylim(lims)
zlim(lims)
view([0 1 0])
camup([1 0 0])
zticks(0)
xticks(0)
yticks(0)
xticklabels([])
zticklabels([])
yticklabels([])
sp2.Box = 'on';

% image scatter plot from positive z-axis: y-axis up
sp3 = subplot(3,1,3);
scatter3(Y2(1,:),Y2(2,:),Y2(3,:),mkrSize,'filled','MarkerEdgeColor',mkrEdgeColor,'MarkerFaceColor',mkrFaceColor,'MarkerFaceAlpha',mkrFaceAlpha,...
    'LineWidth',mkrLineWidth,'MarkerEdgeAlpha',mkrEdgeAlpha)
hold on
grid on
z2 = z2/arrowScale/mx;
scatter3(z2(1),z2(2),z2(3),'filled','Color','r')
xlim(lims)
ylim(lims)
zlim(lims)
view([0 0 1])
camup([0 1 0])
zticks(0)
xticks(0)
yticks(0)
xticklabels([])
zticklabels([])
yticklabels([])
sp3.Box = 'on';

%% VISUALIZE AXIS ESTIMATION FOR METHOD 3

% get gyro data from method 2
Y3 = swing.accel;

% plot only every nth point
i3 = 1:5:size(Y3,2);
Y3 = Y3(:,i3);

% get estimated z axis 
z3 = method3.axis.estimate;

% rotate to body frame
Y3 = R*Y3;
z3 = R*z3;

% set limits for visualization
mx = max(max(abs(Y3)));
lims = limScale*mx*[-1 1];

% stretch z for visualization
z3 = arrowScale*mx*z3;

% init figure
fig2 = figure;
fig2.Position = [1243 83 324 859];

% image scatter plot from positive x-axis: z-axis to left, y-axis up
sp1 = subplot(3,1,1);
scatter3(Y3(1,:),Y3(2,:),Y3(3,:),mkrSize,'filled','MarkerEdgeColor',mkrEdgeColor,'MarkerFaceColor',mkrFaceColor,'MarkerFaceAlpha',mkrFaceAlpha,...
    'LineWidth',mkrLineWidth,'MarkerEdgeAlpha',mkrEdgeAlpha)
grid on
hold on
quiver3(1,0,0,1+z3(1),z3(2),z3(3),'Color',arrowColor,'LineWidth',arrowWidth,'MaxHeadSize',0.3);
xlim(lims)
ylim(lims)
zlim(lims)
view([1 0 0])
camup([0 1 0])
zticks(0)
xticks(0)
yticks(0)
xticklabels([])
zticklabels([])
yticklabels([])
sp1.Box = 'on';

% image scatter plot from positive y-axis: x-axis up
sp2 = subplot(3,1,2);
scatter3(Y3(1,:),Y3(2,:),Y3(3,:),mkrSize,'filled','MarkerEdgeColor',mkrEdgeColor,'MarkerFaceColor',mkrFaceColor,'MarkerFaceAlpha',mkrFaceAlpha,...
    'LineWidth',mkrLineWidth,'MarkerEdgeAlpha',mkrEdgeAlpha)
grid on
hold on
quiver3(0,1,0,z3(1),1+z3(2),z3(3),'Color',arrowColor,'LineWidth',arrowWidth,'MaxHeadSize',0.3);
xlim(lims)
ylim(lims)
zlim(lims)
view([0 1 0])
camup([1 0 0])
zticks(0)
xticks(0)
yticks(0)
xticklabels([])
zticklabels([])
yticklabels([])
sp2.Box = 'on';

% image scatter plot from positive z-axis: y-axis up
sp3 = subplot(3,1,3);
scatter3(Y3(1,:),Y3(2,:),Y3(3,:),mkrSize,'filled','MarkerEdgeColor',mkrEdgeColor,'MarkerFaceColor',mkrFaceColor,'MarkerFaceAlpha',mkrFaceAlpha,...
    'LineWidth',mkrLineWidth,'MarkerEdgeAlpha',mkrEdgeAlpha)
hold on
grid on
z3 = z3/arrowScale/mx;
scatter3(z3(1),z3(2),z3(3),'filled','Color','r')
xlim(lims)
ylim(lims)
zlim(lims)
view([0 0 1])
camup([0 1 0])
zticks(0)
xticks(0)
yticks(0)
xticklabels([])
zticklabels([])
yticklabels([])
sp3.Box = 'on';
