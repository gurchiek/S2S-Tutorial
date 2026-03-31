function [sync,shift,sigrms] = digitalSignalSynchronizationRMS(y1,t1,y2,t2,minShift,maxShift,t1syncStart,t1syncEnd)

% finds time shift of signal 2 necessary to minimize RMS difference between
% signals 1 and 2 given minimum shift minShift and maximum shift maxShift
% (the zero shift is also considered). User can optionally sync based on
% only a portion of signal 1, namely that within the time interval
% t1sycnStart and t2syncStart

% INPUTS
% y1 - 1xN, signal 1
% t1 - 1xN, timestamps for signal 1
% y2 - 1xN, signal 2
% t2 - 1xN, timestamps for signal 2
% minShift - double > 0, minimum non-zero time shift to consider 
% maxShift - double > 0, maximum time shift to consider
% t1syncStart/End - double, see description. Can set to -inf/inf to sync
%                   based on the entire signal 1

% OUTPUTS
% sync.t2 - timestamps for signal 2 synchronized with signal 1
% sync.addThisToSignal2TimeArrayForSynchronization - optimal time shift
% shift - array of timeshifts considered
% sigrms - RMS for each shift considered

% error check
if minShift <= 0
    error('minShift must be non-negative')
end
if maxShift <= minShift
    error('maxShift must be greater than minShift')
end

% get signal 1 within desired time interval
syncIndices = t1 >= t1syncStart & t1 <= t1syncEnd;
t1sync = t1(syncIndices);
y1sync = y1(syncIndices);

% array of shifts to consider
shift = unique(sort([-maxShift:minShift:maxShift,0,maxShift],'ascend'));
nshifts = length(shift);

% for each time shift
sigrms = zeros(1,nshifts);
for k = 1:nshifts

    % shift t2
    tshift = t2 + shift(k);

    % get t1,y1 within tshift
    t1comp = t1sync(t1sync >= tshift(1) & t1sync <= tshift(end));
    y1comp = y1sync(t1sync >= tshift(1) & t1sync <= tshift(end));

    % interpolate at t1 instances within sync window
    y2sync = interp1(tshift,y2,t1comp,'pchip');

    % calc rms
    sigrms(k) = rms(y1comp - y2sync);

end

% get shift that minimizes rms difference
[~,i] = min(sigrms);
sync.t2 = t2 + shift(i);
sync.addThisToSignal2TimeArrayForSynchronization = shift(i);


end