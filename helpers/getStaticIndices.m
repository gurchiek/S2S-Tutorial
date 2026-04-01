function [ indices, v ] = getStaticIndices( data, nsamples, Lnorm)
%Reed Gurchiek, 2017
%   staticndx finds the start and end indices of the interval of length
%   nsamples for the sum of the normed z scores across all data channels is
%   minimum. The norm applied to zscores can be absolute value (norm = 1)
%   or squares (norm = 2)
%
%------------------------------INPUTS--------------------------------------
%
%   data:
%       m x n matrix.  staticndx assumes the largest dimension (m or n)
%       represents indices of observations for different variables, the
%       number of which is equal to the lesser dimension (m or n).  The
%       larger must be greater than nsamples
%
%   nsamples:
%       the number of samples in the static interval desired
%
%   Lnorm:
%       int = 1 or 2, determines whether to sum the absolute values (norm =
%       1) or the squares (norm = 2) of the z-scores
%
%------------------------------OUTPUTS-------------------------------------
%
%   indices:
%       1 x 2.  indices(1) = start index.  indices(2) = end index.
%
%--------------------------------------------------------------------------

%% FIND STATIC

% use abs norm by default
if nargin == 2
    Lnorm = 2;
end

% get dims
[nrows,ncols] = size(data);
if nrows == ncols
    error('data matrix must not be square')
end

% determine dimension of observation indices and number of observations
[N,idim] = max([nrows ncols]);

% verify correct input
if nsamples >= N
    error('nsamples must be larger than largest dimension size of data')
end

% if necessary, transpose data so that each row is a different signal and
% each column corresponds to a discrete time sample
if idim == 1
    data = data';
end

% normalize each signal using z-score based on signal specific mean and
% standard deviation
z = zscore(data,0,2);

% apply user specified norm to the zscores to remove sign (only interested
% in a distance from the mean)
if Lnorm == 1
    z = abs(z);
else
    z = z.^2;
end

% determine number of intervals of length nsamples
nint = N - nsamples + 1;

% allocate space for v where v(k) is a measure of how variable the signal
% was around its own mean within interval k
v = zeros(1,nint);

% for each interval of length nsamples
for k = 1:nint
    
    % first get average normed z score for each signal within interval
    temp = mean(z(:,k:k+nsamples-1),2);

    % now average across all signals
    v(k) = mean(temp);
    
end
    
% get starting index for the most static interval
[~,indices(1)] = min(v);
    
% get ending index for most static interval
indices(2) = indices(1) + nsamples - 1;     

end

