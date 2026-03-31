function shank = reconstructShankMotionFromMarkerData(trcStaticCalibrationData,trcMotionData,originMarkerName)

% solves unconstrained inverse kinematics problem for shank body. The
% origin is taken as the marker specified by originMarkerName (e.g.,
% 'LMAL', which is the default).

%% INIT

% use LMAL as default origin marker
if nargin == 2
    originMarkerName = 'LMAL';
end

% shank marker names
shankMarkers = {'LMAL','FIBHEAD','MMAL','TIBTUB','SHANK1','SHANK2','SHANK3','IMU1','IMU2','MEPI','LEPI'};
nmkr = length(shankMarkers);

%% GET SHANK CONFIGURATION IN STATIC CALIBRATION TRIAL

% loop over all markers in static trial and averag
for k = 1:nmkr
    ref_mkr.(shankMarkers{k}).position = mean(trcStaticCalibrationData.marker.(shankMarkers{k}).position,2);
end

% ISB anatomical locations
LM = ref_mkr.LMAL.position;
MM = ref_mkr.MMAL.position;
IM = (LM + MM) / 2;
MC = ref_mkr.MEPI.position;
LC = ref_mkr.LEPI.position;
IC = (MC + LC) / 2;

% Wu et al 2002 3.2.3
x = normalize(cross(MC - IM, LC - IM),1,'norm'); % frontal
z = normalize(cross(x, IC - IM),1,'norm'); % sagittal
y = normalize(cross(z,x),1,'norm'); % transverse

% orientation in reference config
ref_R = [x y z];
ref_q = convdcm(ref_R,'q');
ref_r = ref_mkr.(originMarkerName).position;

%% GET SHANK CONFIGURATION IN MOTION TRIAL

% get struct of displaced shank marker positions
for k = 1:nmkr
    disp_mkr.(shankMarkers{k}).position = trcMotionData.marker.(shankMarkers{k}).position;
end

% solve unconstrained inverse kinematics
[q,p] = rigidBodyDisplacement_v3(ref_q,ref_r,ref_mkr,disp_mkr);
R = convq(q,'dcm');

% num samples and time array
N = trcMotionData.numSamples;
t = trcMotionData.time;

% time differentiate position
v = fdiff5(p,t);

% time differentiate dcm
Rdot = fdiffmat(R,t,5);

% get angular velocity in body frame
w = zeros(3,N);
for k = 1:N

    % calc angular velocity in skew symmetric form
    wskew = R(:,:,k)' * Rdot(:,:,k);
    wskew = (wskew - wskew') / 2; % ensure skew symmetric

    % unskew
    w(1,k) = wskew(3,2);
    w(2,k) = wskew(1,3);
    w(3,k) = wskew(2,1);

end

% lowpass filter at 6 hz
w = bwfilt(w,6,trcMotionData.samplingFrequency,'low',4);

% store
shank.R = R;
shank.wB = w;
shank.p = p;
shank.v = v;
shank.time = t;
shank.sf = trcMotionData.samplingFrequency;

end