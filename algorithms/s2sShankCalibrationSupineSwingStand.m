function results = s2sShankCalibrationSupineSwingStand(supine,swing,stand)

% INPUTS
% supine, swing, stand - structs with 3 fields: accel, gyro, and sf
%   accel, gyro - 3xn arrays of accel (m/s/s) or gyro (rad/s) data
%   sf - sampling frequency
% supine/stand data are static postures and should satisfy conditions C1
% swing data are uniaxial rotation data and should satisfy conditions C2

% OUTPUTS
% results.stand.calibration - output from estimateAxisFromStaticPosture
% results.supine.calibration - output from estimateAxisFromStaticPosture
% results.swing.calibration - output from estimateAxisFromUniaxialRotation
% results.R_S2B - matrix that rotates from sensor to body frame
% results.gyro - gyro bias and noise variance
% results.accel - accel scale factor, noise variance, bias, and local gravity

% to calibrate accel (ya) and gyro (yw) data use:
% R = results.R_S2B;
% bw = results.gyro.bias;
% c_ = results.accel.inverse_scale_factor
% ba = results.accel.bias
% ya_cal = c_ * R * ya (if you believe null bias)
% ya_cal = R * (ya - ba) (if you believe null scale factor offset)
% yw_cal = R * (yw - bw)

%--------------------------------------------------------------------------
% Copyright 2026 Reed Gurchiek
% 
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
% 
%     http://www.apache.org/licenses/LICENSE-2.0
%--------------------------------------------------------------------------

%% SETTINGS

% local gravitational acceleration
% OpenSim: 9.8066499999999994
% based on clemson elev (221 m) and lat (34.68 deg) = 9.796382780611296
elev = 221; % meters
lat = 34.68; % degrees
g = gravitywgs84(elev,lat);

% options for static posture calibration
options.static_interval_duration = 1; % 1 second of data
options.z_score_magnitude_method = 2;
options.use_gyro_data_for_static_interval_detection = true;
options.g = g;

% axis estimated in standing posture and representation wrt body
results.stand.axis_name = 'y';
yB = [0 1 0]';
results.stand.axis_representation_wrt_body = yB;

% axis estimated in supine posture and representation wrt body
results.supine.axis_name = 'x';
xB = [1 0 0]';
results.supine.axis_representation_wrt_body = xB;

% axis estimated in swing motion and representation wrt body
results.swing.axis_name = 'z';
zB = [0 0 1]';
results.swing.axis_representation_wrt_body = zB;

%% AXIS ESTIMATION

% stand
options.static_interval_samples = round(options.static_interval_duration * stand.sf);
results.stand.calibration = estimateAxisFromStaticPosture(stand,options);
yS = results.stand.calibration.axis.estimate;
wy = 1/results.stand.calibration.axis.uncertainty;

% supine
options.static_interval_samples = round(options.static_interval_duration * supine.sf);
results.supine.calibration = estimateAxisFromStaticPosture(supine,options);
xS = results.supine.calibration.axis.estimate;
wx = 1/results.supine.calibration.axis.uncertainty;

% gyro bias and noise variance in static postures
bw1 = results.stand.calibration.gyro.bias;
bw2 = results.supine.calibration.gyro.bias;
wn1 = results.stand.calibration.gyro.noise_variance;
wn2 = results.supine.calibration.gyro.noise_variance;

% estiamte gyro bias as weighted average, weight by inverse noise variances
results.gyro.bias = (bw1/wn1 + bw2/wn2) / (1/wn1 + 1/wn2); % weighted average based on inverse noise variances
results.gyro.noise_variance = min([wn1 wn2]);

% estimate isotropic scale factor assuming zero accel bias
c1 = results.stand.calibration.accel.nonunity_scale_factor;
c2 = results.supine.calibration.accel.nonunity_scale_factor;
an1 = results.stand.calibration.accel.noise_variance_nonunity_scale;
an2 = results.supine.calibration.accel.noise_variance_nonunity_scale;
c = (c1/an1 + c2/an2) / (1/an1 + 1/an2); % weighted average based on inverse noise variances
results.accel.inverse_scale_factor = 1/c;
results.accel.noise_variance = min([an1 an2]);
results.accel.gravitational_acceleration = g;

% initial guess for z = cross(x,y) using estimated axes from stand/supine
z0 = cross(xS,yS);

% swing
swing.gyro = swing.gyro - results.gyro.bias;
results.swing.calibration = estimateAxisFromUniaxialRotation(swing,struct('initial_guess',z0));
zS = results.swing.calibration.axis.estimate;
wz = 1/results.swing.calibration.axis.uncertainty;

%% ORIENTATION ESTIMATION

% construct matrices with axis representations in body and sensor
Usensor = [xS yS zS];
Ubody = [xB yB zB];

% construct weights vector
w = [wx wy wz];
w = w / sum(w);

% solve wahbas problem
results.wahba = estimateOrientationFromAxisEstimates(Ubody,Usensor,w);
results.R_S2B = results.wahba.R;

%% ESTIMATE OF ACCEL BIAS

% assuming accel scale factor = 1, can estimate accel bias given local
% gravitational acceleration, orientation of posture-specific body-fixed
% axis in sensor frame from results.R_S2B', and assuming that body-fixed
% axis is exactly aligned with gravity

% estimate bias in stand
ybar1 = results.stand.calibration.accel.mean;
uS1 = results.R_S2B' * results.stand.axis_representation_wrt_body;
ba1 = ybar1 - g * uS1;

% estimate bias in supine
ybar2 = results.supine.calibration.accel.mean;
uS2 = results.R_S2B' * results.supine.axis_representation_wrt_body;
ba2 = ybar2 - g * uS2;

% combine in a weighted average
ba = (ba1/an1 + ba2/an2) / (1/an1 + 1/an2);
results.accel.bias = ba;

%% REPORT

fprintf('\n---------------------------------------------------\n')
fprintf('Shank sensor 2 segment calibration results\n\n')

% x-axis
ux = results.R_S2B' * xB;
fprintf('-Estimated shank x-axis in sensor frame: %5.4f, %5.4f, %5.4f\n',ux(1),ux(2),ux(3))
fprintf('     -From supine:  %5.4f, %5.4f, %5.4f (angular error = %4.2f deg)\n',xS(1),xS(2),xS(3),acosd(ux'*xS))
fprintf('           -Weight in Wahba''s problem: %4.3f\n',w(1))

% y-axis
uy = results.R_S2B' * yB;
fprintf('-Estimated shank y-axis in sensor frame: %5.4f, %5.4f, %5.4f\n',uy(1),uy(2),uy(3))
fprintf('     -From stand:  %5.4f, %5.4f, %5.4f (angular error = %4.2f deg)\n',yS(1),yS(2),yS(3),acosd(uy'*yS))
fprintf('           -Weight in Wahba''s problem: %4.3f\n',w(2))

% z-axis
uz = results.R_S2B' * zB;
fprintf('-Estimated shank z-axis in sensor frame: %5.4f, %5.4f, %5.4f\n',uz(1),uz(2),uz(3))
fprintf('     -From swing:  %5.4f, %5.4f, %5.4f (angular error = %4.2f deg)\n',zS(1),zS(2),zS(3),acosd(uz'*zS))
fprintf('           -Weight in Wahba''s problem: %4.3f\n',w(3))

% miscellaneous
fprintf('-Gyro bias estimate: %5.4f, %5.4f, %5.4f\n',results.gyro.bias)
fprintf('     -Magnitude: %5.4f rad/s\n',vecnorm(results.gyro.bias))
fprintf('-Gyro noise variance: %.3e\n',results.gyro.noise_variance)
fprintf('-Accel noise variance: %.3e\n',results.accel.noise_variance)
fprintf('-Accel scale factor (assuming no bias): %4.3f\n',1/results.accel.inverse_scale_factor)
fprintf('-Accel bias estimate (assuming no scale factor offsets): %4.3f, %4.3f, %4.3f\n',results.accel.bias)
fprintf('     -Magnitude: %5.4f m/s/s\n',vecnorm(results.accel.bias))

fprintf('\n---------------------------------------------------\n')

end