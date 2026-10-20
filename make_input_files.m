function make_input_files(cfg)
%MAKE_INPUT_FILES  Create the default input .mat files of the simulator.
%
%   MAKE_INPUT_FILES() writes the four input files named in radar_config.m:
%     antenna_pattern.mat : receive elements, element spacing, azimuth pattern,
%                           elevation element pattern, transmit elevation pattern, gains
%     ddc_filters.mat     : low-pass FIR filters of the two DDC stages (h1, h2)
%     dbf_matrix.mat      : elevation DBF matrix W [nBeams x nElements] and beam elevations
%     stagger_prf.mat     : staggered PRI vector and the matching MTI weights
%
%   Replace any of these files with your own data (same variable names) to
%   simulate a different radar.

if nargin < 1, cfg = radar_config(); end
addpath(fullfile(fileparts(mfilename('fullpath')), 'src'));

lambda = cfg.const.c / cfg.fc;
folder = fileparts(cfg.files.antenna);
if ~exist(folder, 'dir'), mkdir(folder); end

%% 1) Antenna
nElements        = 48;          % receive elements (rows) in elevation
elementSpacing_m = lambda / 2;  % element spacing in elevation [m]
azBeamwidth_deg  = 4.1;         % 3 dB azimuth beamwidth (one way)
elBeamwidth_deg  = 70;          % 3 dB elevation beamwidth of a single element (row)

% Azimuth pattern: Taylor-weighted horizontal line array scaled to the requested beamwidth
az_deg   = (-180:0.02:180)';
nAz      = 32;
wAz      = taylor_window(nAz, 4, -30);
afPower  = @(d, az) abs(exp(1j*2*pi*d*(0:nAz-1)'*sind(az(:)')).' * wAz).^2 / sum(wAz)^2;
halfBw   = @(d) 10*log10(afPower(d, azBeamwidth_deg/2)) + 3.0103;   % zero at the right spacing
dAz      = fzero(halfBw, [0.3 0.7]);                                % element spacing [wavelengths]
azGain_dB = max(10*log10(afPower(dAz, az_deg)), -50);
azGain_dB(abs(az_deg) > 90) = -50;                                  % no back lobe

% Elevation element pattern cos^q(theta), theta measured from the array broadside
elRel_deg     = (-90:0.1:90)';
q             = log(0.5) / log(cosd(elBeamwidth_deg/2));             % 3 dB at +-35 deg
elElemGain_dB = max(10*log10(cosd(elRel_deg).^q + eps), -30);
txElGain_dB   = elElemGain_dB;   % the transmitter uses the same wide elevation beam

% Peak gains (rule of thumb G = 26000 / (az_bw * el_bw) with beamwidths in degrees)
rxChannelGain_dBi = 10*log10(26000 / (azBeamwidth_deg * elBeamwidth_deg));
txGain_dBi        = rxChannelGain_dBi;

save(cfg.files.antenna, 'nElements', 'elementSpacing_m', 'az_deg', 'azGain_dB', ...
    'elRel_deg', 'elElemGain_dB', 'txElGain_dB', 'rxChannelGain_dBi', 'txGain_dBi', '-v7');

%% 2) DDC low-pass filters (windowed sinc, Blackman window, unit DC gain)
fs  = cfg.adc.fs_Hz;
fs1 = fs / cfg.ddc.decimation(1);
h1  = lowpass_fir(70, 5.8e6 / fs);    % stage 1: pass +-1 MHz, stop above 10.5 MHz (fs = 120 MHz)
h2  = lowpass_fir(66, 1.0e6 / fs1);   % stage 2: pass +-0.5 MHz, stop above 1.5 MHz (fs = 12 MHz)
save(cfg.files.ddc, 'h1', 'h2', '-v7');

%% 3) Elevation DBF matrix: 21 beams equally spaced in u = sin(angle from broadside)
nBeams     = 21;
uFirst     = sind(0.5 - cfg.ant.tilt_deg);       % lowest beam 0.5 deg above the horizon
uStep      = 0.04;                               % about one beamwidth
beamU      = uFirst + (0:nBeams-1)' * uStep;
beamEl_deg = cfg.ant.tilt_deg + asind(beamU);    % beam elevations (at the azimuth boresight)
taper      = taylor_window(nElements, 4, -30);
n          = 0:nElements-1;
W          = (taper(:).' .* exp(-1j*2*pi*elementSpacing_m/lambda * beamU * n));
W          = W ./ sqrt(sum(abs(W).^2, 2));       % unit noise gain for every beam
save(cfg.files.dbf, 'W', 'beamEl_deg', '-v7');

%% 4) Staggered PRF: PRI vector and time-varying 3-pulse MTI weights
priVector_s = [920; 1000; 960; 1040] * 1e-6;
K = numel(priVector_s);
mtiCoeffs = zeros(K, 3);
for k = 1:K
    % Pulses k, k+1, k+2 are separated by T1 and T2. The weights null DC and
    % the first-order (linear phase) clutter term; they equal [1 -2 1] when T1 = T2.
    T1 = priVector_s(k);
    T2 = priVector_s(mod(k, K) + 1);
    mtiCoeffs(k, :) = 2 * [T2, -(T1 + T2), T1] / (T1 + T2);
end
save(cfg.files.stagger, 'priVector_s', 'mtiCoeffs', '-v7');

fprintf('Input files written to %s\n', folder);
end

function h = lowpass_fir(nTaps, fCutNorm)
% Windowed-sinc low-pass FIR. fCutNorm = cut-off frequency / sampling frequency.
n = (0:nTaps-1)' - (nTaps-1)/2;
x = 2*fCutNorm*n;
s = ones(size(x));
s(x ~= 0) = sin(pi*x(x ~= 0)) ./ (pi*x(x ~= 0));
k = (0:nTaps-1)';
blackman = 0.42 - 0.5*cos(2*pi*k/(nTaps-1)) + 0.08*cos(4*pi*k/(nTaps-1));
h = s .* blackman;
h = h / sum(h);
end
