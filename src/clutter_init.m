function clut = clutter_init(sys)
%CLUTTER_INIT  Static ground clutter map around the radar.
%   clut = CLUTTER_INIT(sys) builds a grid of clutter cells (azimuth x range)
%   with complex Gaussian (Rayleigh amplitude) reflectivity of constant sigma0,
%   a random internal-motion Doppler per cell, and precomputes what
%   clutter_signal.m needs to synthesise the clutter echo of every PRI.

cfg = sys.cfg;
cc  = cfg.clutter;
c   = cfg.const.c;
fs  = sys.fs;

% Range grid: cell size is a whole number of ADC samples so delays fall on samples
dR = cc.rangeStepSamples * c / (2*fs);
rMax = cc.maxRange_m;
if isempty(rMax)
    rMax = sqrt(2 * 4/3 * 6371e3 * cfg.ant.height_m);       % radar horizon (4/3 earth)
end
clut.range_m = (ceil(cc.minRange_m/dR):floor(rMax/dR)) * dR;
clut.delaySamples = round(2*clut.range_m/c*fs);           % two-way delay [ADC samples]
clut.az_deg = (0:cc.azStep_deg:360-cc.azStep_deg)';
nAz = numel(clut.az_deg);
nR  = numel(clut.range_m);

% Mean power of a cell at the ADC (antenna boresight in azimuth) and its rms IF amplitude
[~, u] = steering_vector(sys, 0, cc.elevation_deg);
elRel  = asind(u);
gEl    = pattern_gain(sys.ant.elRel_deg, sys.ant.txElGain_dB, elRel) ...
       * pattern_gain(sys.ant.elRel_deg, sys.ant.elElemGain_dB, elRel);
area   = clut.range_m * deg2rad(cc.azStep_deg) * dR;      % cell area [m^2]
power  = sys.radarConst * 10^(cc.sigma0_dB/10) * area ./ clut.range_m.^4 * gEl;  % [W]
ampRms = sqrt(2*power);                                   % rms IF amplitude at the receiver input [sqrt(W)]
clut.rho = ampRms .* (randn(nAz, nR) + 1j*randn(nAz, nR)) / sqrt(2);
clut.fd_Hz = 2*cc.icmStd_mps/sys.lambda * randn(nAz, nR); % internal clutter motion

% IF transmit waveforms at the ADC rate and their FFTs (for the fast convolution)
names = {'short', 'long'};
for i = 1:2
    wv = sys.wave.(names{i});
    nTx = round(wv.width_s * fs);
    t   = (0:nTx-1)' / fs;
    tx  = waveform_envelope(wv, t) .* exp(1j*2*pi*sys.fIF*t);
    nFft = 2^nextpow2(sys.nAdc + nTx - 1);
    clut.(names{i}).nTx    = nTx;
    clut.(names{i}).nFft   = nFft;
    clut.(names{i}).txFft  = fft(tx, nFft);
    % ADC sample (0-based, relative to the window start) where an echo of zero delay starts
    clut.(names{i}).start0 = round((wv.tTx_s - sys.rxStart_s) * fs);
end
end
