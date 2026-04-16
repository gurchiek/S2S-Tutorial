function cal = estimateAxisFromPlanarMotion(data,options)

% estimate rotation axis in sensor frame from accelerometer data assuming 
% planar motion conditions C3

% INPUTS
% data.accel - 3xN array of accelerometeter data in m/s/s
% options.initial_guess - 3x1 vector, initial guess of the estimated axis
% options.store_calibration_data - boolean, default = false

% OUTPUTs
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

%% ALGORITHM 3

% unpack
Ya = data.accel;

% options
if nargin == 1
    options = struct;
end
if ~isfield(options,'store_calibration_data')
    options.store_calibration_data = false;
end

% num samples
N = size(Ya,2);

% calculate sample average
ybar = sum(Ya,2) / N;

% calculate scaled covariance matrix
C = (Ya * Ya')/N - ybar * ybar';

% eigen decomp
[V,lam] = eig(C,'nobalance','vector');
[lam,isort] = sort(lam,'ascend');
V = V(:,isort);

% estimate rotation axis
u = V(:,1);

% correct for sign if initial guess available
if isfield(options,'initial_guess')
    u = sign(u' * options.initial_guess) * u;
end

% estimate mean dot product between u and accel measurements = projection
% of gravity onto u + bias
mu = ybar' * u;

% noise variance estimate
Ess = N * lam(1);
sigsqr = lam(1);

% can verify Ess = lam(1) is equal to the following
% sum((u'*Ya-u'*ybar).^2) 
% u'* N*C * u 

% axis estimate uncertainty
uvar = sigsqr * (1/lam(2) + 1/lam(3)) / N;

% can verify uvar with following:
% Iu = N*C / sigma2; % unconstrained observation Fisher info matrix
% Z = V(:,2:end); % basis of null space of constraint jac
% CRLB = Z / (Z' * Iu * Z) * Z';
% uvar = trace(CRLB);

%% PACKUP OUTPUTS

% axis estimate
cal.axis.estimate = u;
cal.axis.uncertainty = uvar;
cal.axis.sign_determinant = isfield(options,'initial_guess');

% mu = g * u' * yW_S = projection of gravity onto = constant under conditions C3 up to noise
cal.mu = mu; 

% stats
cal.stats.num_samples = N;
cal.stats.sse = Ess;
cal.stats.mean_angular_error_deg = mean(asind(abs(u' * normalize(Ya - mu,1,'norm'))));

% aceelerometer-specific
if options.store_calibration_data
    cal.accel.data = Ya;
end
cal.accel.noise_variance = sigsqr;
cal.accel.mean = ybar;
cal.accel.covariance.matrix = C;
cal.accel.covariance.evals = lam;
cal.accel.covariance.evecs = V;

end