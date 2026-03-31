function stride = detectGaitEventsShankGyro(t,w,sf)

% event detection algo from Mansour et al. 2015: swing  phase associated
% with large positive peaks in angular rate signal. These are identified
% first after low pass filtering at 3 hz. 
% Before each peak is a zero crossing and the foot off event is
% associated with the negative peak in the raw gyro signal just before this
% crossing. Likewise, after each peak is a zero crossing and the foot
% contact event is associaged with the negative peak in the raw gyro signal
% just after this crossing. 

% gait classification not perfect and neither is event detection, some
% non-stride data may make it through. Best practice is to apply
% lower/upper bounds on gait params (e.g. stride time, duty factor, etc)
% and remove outliers. This done in walkingStrideAnalysis_basic

% INPUTS:
% t - 1xn, timestamps
% w - 3xn, gyro (lateral distal shank), rad/s
% a - 3xn, accel (lateral distal shank), g

% OUTPUTS
% stride, 1xm struct, fields specify:
%   (1) timestamp of swing phase
%   (2) timestamp of previous foot off event
%   (3) timestamp of next foot contact event
%   (4) timestamp of the start of the walking bout this stride is in (t(1))
%   (5) timestamp of the end of the walking bout this stride is in (t(end))

warning('off')

%% parameter settings

swing_gyro_min_peak = 1.5;

%% get gait events

% handle case that user input the z axis only already. Code operates by
% default on 3xn gyro data but needs only the z axis (w(3,:)). If user gave
% w(3,:) as the input then simply make a 3xn array with the first and
% second rows all zero
[nrows,ncols] = size(w);
if nrows > ncols
    w = w'; 
    [nrows,ncols] = size(w);
end
if nrows == 1
    w = [zeros(2,ncols); w];
end


% low pass gyro @ 3
wz = bwfilt(w(3,:),3,sf,'low',4);
nsamp = length(wz);

% get all positive peaks > swingpeakmin rad/s
[~,ipos] = findpeaks(wz,'MinPeakHeight',swing_gyro_min_peak);
nEstimatedStrides = length(ipos);

% get all negative peaks on raw sagittal plane signal
[~,ineg] = findpeaks(-w(3,:),'MinPeakHeight',0);

% get all positive and negative going zero crossings
[i,type] = getZeroCrossings(wz);
pz0 = i(contains(type,{'n2p','n2z','z2p'}));
nz0 = i(contains(type,{'p2n','p2z','z2n'}));

% loop over each positive peak (ipos), there should be a positive going 0
% crossing just before this and a negative going one just after this, these
% are the only ones we want to keep for event detection, others could occur
% if the angular rate signal goes positive during stance, which is rare but
% happens on occasion
pz = zeros(1,nEstimatedStrides);
nz = zeros(1,nEstimatedStrides);
iremovepz = false(1,nEstimatedStrides);
iremovenz = false(1,nEstimatedStrides);
for k = 1:nEstimatedStrides
    
    % get all pos going crossing prior to this mid swing peak
    x = pz0(pz0 < ipos(k));

    % keep only the last one
    if ~isempty(x)
        pz(k) = x(end);
    else
        iremovepz(k) = true;
    end

    % get all neg going crossings after this mid swing peak
    x = nz0(nz0 > ipos(k));

    % keep only the first one
    if ~isempty(x)
        nz(k) = x(1);
    else
        iremovenz(k) = true;
    end
    
end
pz(iremovepz) = [];
nz(iremovenz) = [];

% initialize stride struct
stride(1:nEstimatedStrides) = struct('swingPhaseTimestamp',[],...
                                     'firstFootOffTimestamp',[],...
                                     'firstNegativeGoing0CrossingTimestamp',[],...
                                     'firstFootContactTimestamp',[],...
                                     'secondNegativeGoing0CrossingTimestamp',[],...
                                     'secondFootContactTimestamp',[],...
                                     'secondPositiveGoing0CrossingTimestamp',[],...
                                     'secondFootOffTimestamp',[],...
                                     'boutStartTimestamp',t(1),...
                                     'boutEndTimestamp',t(end));

% for each positive peak (swing phase)
for k = 1:nEstimatedStrides
    
    % timestamp associated with peak angular velocity in swing
    stride(k).swingPhaseTimestamp = t(ipos(k));
    
    % get positive going crossings before mid swing
    x = pz(pz < ipos(k));
    if ~isempty(x)

        % interpolate and store last one
        % stride(k).firstPositiveGoing0CrossingTimestamp = t(x(end)); % non-interpolated answer
        if x(end)-2 >= 1
            interpInd = x(end)-2:x(end)+1;
            interpMethod = 'pchip';
        else
            interpInd = x(end)-1:x(end)+2;
            interpMethod = 'spline';
        end
        stride(k).firstPositiveGoing0CrossingTimestamp = interp1(wz(interpInd),t(interpInd),0,interpMethod);

        % get negative peaks in raw signal before this
        x = ineg(ineg < x(end));

        % store last one as foot off
        if ~isempty(x); stride(k).firstFootOffTimestamp = t(x(end)); end

    end

    % get negative going crossings before mid swing
    x = nz(nz < ipos(k));
    if ~isempty(x)

        % interpolate and store last one
        % stride(k).firstNegativeGoing0CrossingTimestamp = t(x(end)); % non-interpolated answer
        if x(end)-2 >= 1
            interpInd = x(end)-2:x(end)+1;
            interpMethod = 'pchip';
        else
            interpInd = x(end)-1:x(end)+2;
            interpMethod = 'spline';
        end
        stride(k).firstNegativeGoing0CrossingTimestamp = interp1(wz(interpInd),t(interpInd),0,interpMethod);

        % get negative peaks in raw signal after this
        x = ineg(ineg > x(end));

        % store first one as foot contact (stride start)
        if ~isempty(x); stride(k).firstFootContactTimestamp = t(x(1)); end

    end
    
    % get negative going crossings after mid swing
    x = nz(nz > ipos(k));
    
    % continue if non-empty
    if ~isempty(x)

        % interpolate and store first one
        % stride(k).secondNegativeGoing0CrossingTimestamp = t(x(1)); % non-interpolated answer
        if x(1)+1 <= nsamp
            interpInd = x(1)-2:x(1)+1;
            interpMethod = 'pchip';
        else
            interpInd = x(1)-3:x(1);
            interpMethod = 'spline';
        end
        stride(k).secondNegativeGoing0CrossingTimestamp = interp1(wz(interpInd),t(interpInd),0,interpMethod);
        
        % get negative peaks in raw signal after this
        x = ineg(ineg > x(1));

        % store first one as next foot contact timestamp (stride end)
        if ~isempty(x); stride(k).secondFootContactTimestamp = t(x(1)); end
            
    end

    % get positive going zero crossings after mid swing
    x = pz(pz > ipos(k));
    if ~isempty(x)

        % interpolate and store first one
        % stride(k).secondPositiveGoing0CrossingTimestamp = t(x(1)); % non-interpolated answer
        if x(1)+1 <= nsamp
            interpInd = x(1)-2:x(1)+1;
            interpMethod = 'pchip';
        else
            interpInd = x(1)-3:x(1);
            interpMethod = 'spline';
        end
        stride(k).secondPositiveGoing0CrossingTimestamp = interp1(wz(interpInd),t(interpInd),0,interpMethod);

        % get negative peaks in raw signal before this
        x = ineg(ineg < x(1));

        % store last one as next foot off event
        if ~isempty(x); stride(k).secondFootOffTimestamp = t(x(end)); end

    end
    
end

% visualize
% figure
% hold on
% plot(t,wz,'k')
% plot(t,w(3,:),'Color',[0.5 0.5 0.5])
% for k = 1:length(stride)
%     sw = stride(k).swingPhaseTimestamp;
%     scatter(sw,wz(t == sw),'r');
%     fc = stride(k).secondFootContactTimestamp;
%     if ~isempty(fc)
%         scatter(fc,w(3,t == fc),'b')
%     end
%     fo = stride(k).firstFootOffTimestamp;
%     if ~isempty(fo)
%         scatter(fo,w(3,t == fo),'b')
%     end
%     nz1 = stride(k).firstNegativeGoing0CrossingTimestamp;
%     if ~isempty(nz1)
%         scatter(nz1,zeros(1,length(nz1)),'g')
%     end
%     nz2 = stride(k).secondNegativeGoing0CrossingTimestamp;
%     if ~isempty(nz2)
%         scatter(nz2,zeros(1,length(nz2)),'g')
%     end
%     pz1 = stride(k).firstPositiveGoing0CrossingTimestamp;
%     if ~isempty(pz1)
%         scatter(pz1,zeros(1,length(pz1)),'g')
%     end
%     pz2 = stride(k).secondPositiveGoing0CrossingTimestamp;
%     if ~isempty(pz2)
%         scatter(pz2,zeros(1,length(pz2)),'g')
%     end
% end
    
end