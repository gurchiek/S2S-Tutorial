function [q, t] = qinterp(q0,t0,t)

% INPUTS
% q0 - 4xn array of quaternions
% t0 - 1xn time array, should be strictly ascending
% t - 1xm new time array, should be strictly ascending

% OUTPUTS
% q - 4xm array of interpolated quaternions
% t - 1xm new time array, same as input t unless t was input not non
% ascending order and/or non unique

%% qinterp

% make sure t0 and t are unique and ascending
t = sort(unique(t),'ascend');
[t0,i] = unique(t0);
q0 = q0(:,i);
[t0,i] = sort(t0,'ascend');
q0 = q0(:,i);

% throw error if extrapolating
if min(t) < t0(1) || max(t) > t0(end)
    error('Cannot extrapolate: at least one interpolation grid point lies outside the input time grid')
end

% for each new time
q = zeros(4,length(t));
for k = 1:length(t)
    
    % get indices in t0 that bound tk
    i1 = find(t(k) >= t0(1:end-1));
    i1 = i1(end);
    i2 = i1+1;
    
    % get quaternions at interval bounds
    q1 = q0(:,i1);
    q2 = q0(:,i2);
    
    % normalized time step
    tnorm = (t(k) - t0(i1)) / (t0(i2) - t0(i1));
    
    % transition quaternion from q1 to q2
    theta0 = acos(q1'*q2);
    theta = tnorm * theta0;
    
    % interpolate
    if abs(theta0) < eps
        q(:,k) = q1;
    else
        q(:,k) = (sin(theta0-theta)*q1 + sin(theta)*q2) / sin(theta0);
    end
    
end

end