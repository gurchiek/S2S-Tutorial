function cal = estimateAxisFromStaticPosture(data,options)

% estimate direction of gravity in sensor frame from accelerometer data
% assuming static posture conditions C1

% INPUTS
% data.accel = 3xN array of accel data in m/s/s
% data.gyro = (optional) 3xn gyro data in rad/s
% options.static_interval_samples = int, default = inf
% options.use_gyro_data_for_static_interval_detection = bool, default = false 
% options.z_score_magnitude_method = int, 1 or 2 (default), see getStaticIndices
% options.store_calibration_data = bool, default = false
% options.g = local gravitational acceleration, default = 9.8066499999999994

% OUTPUTS
% see PACKUP OUTPUTS section below

%--------------------------------------------------------------------------
% Copyright 2026 Reed Gurchiek
% 
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
% 
%     http://www.apache.org/licenses/LICENSE-2.0
%--------------------------------------------------------------------------

%% ALGORITHM 1

% default options
if nargin == 1
    options = struct;
end
if ~isfield(options,'static_interval_samples')
    options.static_interval_samples = inf;
end
if ~isfield(options,'use_gyro_data_for_static_interval_detection')
    options.use_gyro_data_for_static_interval_detection = false;
end
if ~isfield(options,'z_score_magnitude_method')
    options.z_score_magnitude_method = 2;
end
if ~isfield(options,'store_calibration_data')
    options.store_calibration_data = false;
end
if ~isfield(options,'g')
    options.g = 9.8066499999999994;
end
if ~isfield(data,'gyro')
    data.gyro = [];
end

% unpack
Ya = data.accel;
Yw = data.gyro;
use_gyro = options.use_gyro_data_for_static_interval_detection;
g = options.g;

% gyro data available? if so, then will compute some gyro relevant stats
has_gyro = ~isempty(Yw);

% get interval within which motion was most static
Nstatic = options.static_interval_samples; % num samples requested for calibration
Ndata = size(Ya,2); % num samples in dataset
if Nstatic < Ndata
    N = Nstatic;
    if use_gyro
        ind = getStaticIndices([Ya; Yw], Nstatic, 2);
    else
        ind = getStaticIndices(Ya, Nstatic, 2);
    end
else
    N = Ndata;
    ind = [1 N];
end
ind = ind(1):ind(2);

% get accelerometer and gyroscope data within static interval
Ya = Ya(:,ind);
if has_gyro
    Yw = Yw(:,ind);
end

% get average acceleration
ybar = mean(Ya,2);

% ML estimate of g
xg = sqrt(ybar' * ybar);

% ML estimate of u
u = ybar / xg;

% ML estimate of scale factor c assuming a non-unity isotropic scale factor
c = xg / g;

% noise variance estimate assuming zero bias and unity scale factor 
err1 = Ya - g*u;
Ess1 = sum(dot(err1,err1));
sigsqr1 = Ess1 / 3 / N;

% noise variance estimate assuming zero bias and non-unity isotropic scale
% factor
err2 = Ya - xg*u;
Ess2 = sum(dot(err2,err2));
sigsqr2 = Ess2 / 3 / N;

% axis estimate uncertainty assuming zero bias and unity scale factor
uncertainty1 = 2 * sigsqr1 / (g*g) / N;

% axis estimate uncertainty assuming zero bias and non-unity isotropic
% scale factor
uncertainty2 = 2 * sigsqr2 / (xg*xg) / N;

% other statistics
avar = var(Ya,0,2); % axis specific sample variance
if has_gyro
    wvar = var(Yw,0,2);
    bw = mean(Yw,2); % bias
    werr = Yw - bw;
    wEss = sum(dot(werr,werr));
    wsigsqr = wEss / 3 / N;
end

%% PACKUP OUTPUTS

% axis estimate
cal.axis.estimate = u;
cal.axis.uncertainty = uncertainty1;
cal.axis.uncertainty_nonunity_scale = uncertainty2;

% gravity
cal.g.true = g;
cal.g.estimate = xg;

% stats
cal.stats.num_samples = N;
cal.stats.sse = Ess1;
cal.stats.sse_nonunity_scale = Ess2;
cal.stats.mean_angular_error_deg = mean(acosd(u' * normalize(Ya,1,'norm')));

% accelerometer specific
if options.store_calibration_data
    cal.accel.calibration_data = Ya;
end
cal.accel.noise_variance = sigsqr1;
cal.accel.noise_variance_nonunity_scale = sigsqr2;
cal.accel.nonunity_scale_factor = c;
cal.accel.axis_specific_variance = avar;
cal.accel.mean = ybar;

% gyro specific 
if has_gyro
    if options.store_calibration_data
        cal.gyro.calibration_data = Yw;
    end
    cal.gyro.noise_variance = wsigsqr;
    cal.gyro.axis_specific_variance = wvar;
    cal.gyro.bias = bw;
end

end