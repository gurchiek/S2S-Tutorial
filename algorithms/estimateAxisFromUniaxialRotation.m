function cal = estimateAxisFromUniaxialRotation(data,options)

% estimate axis of rotation represented in sensor frame from gyroscope data
% assuming uniaxial rotation conditions C2

% INPUTS
% data.gyro - 3xN array of gyroscope data in rad/s
% options.initial_guess - 3x1 vector, initial guess of the estimated axis
% options.store_calibration_data = bool, default = false

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

%% ALGORITHM 2

% unpack
Yw = data.gyro;

% options
if nargin == 1
    options = struct;
end
if ~isfield(options,'store_calibration_data')
    options.store_calibration_data = false;
end

% num samples
N = size(Yw,2);

% calculate scatter matrix
A = Yw * Yw';

% eigen decomp
[V,lam] = eig(A,'nobalance','vector');
[lam,isort] = sort(lam,'ascend');
V = V(:,isort);

% estimate rotation axis
u = V(:,3);

% correct for sign if initial guess available
if isfield(options,'initial_guess')
    u = sign(u' * options.initial_guess) * u;
end

% estimate speeds
s = u' * Yw; % note s'*s = lam(3)

% noise variance estimate
Ess = lam(1) + lam(2);
sigsqr = Ess / 3 / N;

% can verify Ess calculation with following:
% err = zeros(3,N);
% for k = 1:N
%     err(:,k) = Yw(:,k) - s(k) * u;
% end
% sse = sum(dot(err,err));

% axis estimate uncertainty
uvar = 2 * sigsqr / lam(3);

% can verify uvar with following:
% Iu = lam(3) / sigma2; % unconstrained observation Fisher info matrix, also = s'*s / sigma2
% Z = V(:,1:2); % basis of null space of constraint jac
% CRLB = Z * inv(Z' * Iu * Z) * Z';
% uvar = trace(CRLB);

%% PACKUP OUTPUTS

% axis estimate
cal.axis.estimate = u;
cal.axis.uncertainty = uvar;
cal.axis.sign_determinant = isfield(options,'initial_guess');

% speeds
cal.speeds = s;

% stats
cal.stats.num_samples = N;
cal.stats.sse = Ess;
cal.stats.mean_angular_error_deg = mean(acosd(abs(u' * normalize(Yw,1,'norm'))));

% gyro specific
if options.store_calibration_data
    cal.gyro.data = Yw;
end
cal.gyro.noise_variance = sigsqr;
cal.gyro.mean = mean(Yw,2);
cal.gyro.scatter.matrix = A;
cal.gyro.scatter.evals = lam;
cal.gyro.scatter.evecs = V;

end