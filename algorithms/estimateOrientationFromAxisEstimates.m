function cal = estimateOrientationFromAxisEstimates(UA,UB,w)

% estimate rotation matrix R that minimizes sum of squared norm of errors
% UA(:,k) - R*UB(:,k) weighted by w(k). Solves using SVD solution to
% Wahba's problem

% INPUTS
% UA - 3xN matrix of vectors represented in frame A
% UB - 3xN matrix of vectors represented in frame B
% w - 1xN array of weights

% OUTPUTS
% see PACKUP OUTPUTS section below
% estimated orientation = cal.R rotates vectors represented in frame B to
% their representation in frame A: vA = cal.R * vB

%--------------------------------------------------------------------------
% Copyright 2026 Reed Gurchiek
% 
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
% 
%     http://www.apache.org/licenses/LICENSE-2.0
%--------------------------------------------------------------------------

%% ALGORITHM 4

% num samples
N = size(UA,2);

% ensure columns of Ux and UB are all unit length
UA = normalize(UA,1,'norm');
UB = normalize(UB,1,'norm');

% normalize weights to ensure unity gain
if nargin == 3
    w = w / sum(w);
else
    w = ones(1,N);
end

% form the weighted attitude matrix
M = zeros(3);
for k = 1:N
    M = M + UA(:,k) * UB(:,k)' * w(k);
end

% svd
[VL,~,VR] = svd(M);

% handedness correction
h = det(VL) * det(VR);
D = diag([1,1,h]);

% estimate R
R = VL * D * VR';

% error stats
err = UA - R * UB;
err2 = dot(err,err);
weighted_sse = dot(w,err2);
unweighted_sse = sum(err2);
angular_errors = acosd(dot(UA,R*UB));

%% PACKUP OUTPUTS

% rotation matrix such that vA = R * vB
cal.R = R;

% stats
cal.stats.num_axes = N;
cal.stats.sse.weighted = weighted_sse;
cal.stats.sse.unweighted = unweighted_sse;
cal.stats.angular_errors_deg = angular_errors;
cal.stats.mean_angular_error_deg = mean(angular_errors);

end
