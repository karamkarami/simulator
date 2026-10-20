function sys = setup_radar(cfg)
%SETUP_RADAR  Load the input files and derive every dependent radar quantity.
%
%   sys = SETUP_RADAR(cfg) returns a struct with the antenna, noise, timing,
%   waveform, filter, axis, CPI and CFAR parameters used by the processing
%   blocks. The raw configuration is kept in sys.cfg.

c  = cfg.const.c;
fs = cfg.adc.fs_Hz;
sys.cfg    = cfg;
sys.lambda = c / cfg.fc;
sys.fs     = fs;
sys.fIF    = cfg.adc.fIF_Hz;                                  % IF carrier at the ADC [Hz]
sys.fLO    = cfg.fc - sys.fIF;                                % receiver local oscillator [Hz]

%% Input files (defaults are built when missing or when made for other cfg values)
ant  = load_input(cfg, 'antenna', struct('fc', cfg.fc));
ddcF = load_input(cfg, 'ddc', struct('fs_Hz', fs, 'decimation', cfg.ddc.decimation));

%% Antenna
sys.ant = ant;
sys.ant.height_m = cfg.ant.height_m;
sys.ant.tilt_deg = cfg.ant.tilt_deg;
sys.ant.rate_dps = cfg.ant.rpm * 360 / 60;                    % rotation rate [deg/s]
sys.ant.azBeamwidth_deg = beamwidth_3dB(ant.az_deg, ant.azGain_dB);
sys.ant.azSpan_deg = 3 * sys.ant.azBeamwidth_deg;             % echoes outside +-azSpan are ignored
sys.nCh = cfg.ant.nElements;

% Received power per element = radarConst * rcs / R^4 * (pattern factors)
sys.radarConst = cfg.tx.peakPower_W * 10^(ant.txGain_dBi/10) * 10^(ant.rxChannelGain_dBi/10) ...
               * sys.lambda^2 / ((4*pi)^3 * 10^(cfg.rx.losses_dB/10));

%% Noise level and ADC scaling (signals are simulated in ADC units = LSB)
sys.adc.fullScale = 2^(cfg.adc.nBits - 1);
sinad      = 10^((6.02*cfg.adc.enob + 1.76)/10);              % full-scale sine / ADC noise
adcNoise2  = sys.adc.fullScale^2/2 / sinad;                   % total ADC noise power [LSB^2]
sys.adc.sigmaAdc = sqrt(max(adcNoise2 - 1/12, 0));            % ADC noise added before rounding [LSB]
sys.adc.sigmaTh  = sqrt(adcNoise2 * 10^(cfg.adc.noiseAboveAdc_dB/10));   % thermal noise [LSB]
sys.N0    = cfg.const.kB * cfg.const.T0 * 10^(cfg.rx.noiseFigure_dB/10); % noise density [W/Hz]
sys.N0eff = sys.N0 * (sys.adc.sigmaTh^2 + adcNoise2) / sys.adc.sigmaTh^2; % incl. ADC noise
sys.adc.scale = sys.adc.sigmaTh / sqrt(sys.N0 * fs/2);        % amplitude scale [LSB / sqrt(W)]

%% PRI sequence and MTI weights
switch lower(cfg.prf.mode)
    case 'fixed'
        sys.pri = cfg.prf.pri_s;
        sys.mtiCoeffs = [1 -2 1];                             % three-pulse canceller
    case 'staggered'
        st = load_input(cfg, 'stagger', []);
        sys.pri = st.priVector_s(:);
        sys.mtiCoeffs = st.mtiCoeffs;                         % one row per stagger position
    otherwise
        error('Unknown PRF mode "%s"', cfg.prf.mode);
end
sys.staggered = strcmpi(cfg.prf.mode, 'staggered');
sys.nStagger  = numel(sys.pri);
if size(sys.mtiCoeffs, 1) == 1
    sys.mtiCoeffs = repmat(sys.mtiCoeffs, sys.nStagger, 1);
end

%% DDC filters and the numerically controlled oscillator (NCO)
D1 = cfg.ddc.decimation(1);
D2 = cfg.ddc.decimation(2);
sys.ddc.h1 = ddcF.h1(:);
sys.ddc.h2 = ddcF.h2(:);
sys.ddc.decim = [D1 D2];
L1 = numel(sys.ddc.h1);
L2 = numel(sys.ddc.h2);
fs1 = fs / D1;
sys.fsOut = fs1 / D2;
% Single filter at the ADC rate equivalent to the two stages (for the noise bandwidth)
h2Up = zeros(D1*(L2-1)+1, 1);
h2Up(1:D1:end) = sys.ddc.h2;
sys.ddc.hEq = conv(sys.ddc.h1, h2Up);
sys.ddc.noiseBw_Hz = fs * sum(sys.ddc.hEq.^2);                % complex noise bandwidth [Hz]

%% PRI timeline: short pulse, receive, long pulse, receive
% The receiver is blanked while a pulse is transmitted. One ADC record per
% PRI runs from the end of the short pulse to the end of the long-pulse
% receive window (it contains the blanked interval of the long pulse).
names = {'short', 'long'};
for i = 1:2
    wv = cfg.wave.(names{i});
    wv.name     = names{i};
    wv.lambda_m = sys.lambda;                                 % both sub-pulses are on fc
    sys.wave.(names{i}) = wv;
end
ws = sys.wave.short;
wl = sys.wave.long;
if wl.tTx_s < ws.tTx_s + ws.width_s
    error('The long pulse must start after the end of the short pulse.');
end
filterSpan    = (L1-1)/fs + (L2-1)/fs1;                       % input span of one DDC output sample
sys.rxStart_s = ceil((ws.tTx_s + ws.width_s)*fs) / fs;        % end of the short pulse
rxEnd = wl.tTx_s + wl.width_s + 2*cfg.rx.maxRange_m/c + filterSpan + 2/sys.fsOut;
if rxEnd > min(sys.pri)
    error('Receive window (%.1f us) does not fit in the shortest PRI; reduce cfg.rx.maxRange_m.', rxEnd*1e6);
end
sys.nAdc = round((rxEnd - sys.rxStart_s) * fs);
sys.tAdc = sys.rxStart_s + (0:sys.nAdc-1)' / fs;             % ADC sample times in the PRI [s]
sys.rxGate = single(sys.tAdc < wl.tTx_s | sys.tAdc >= wl.tTx_s + wl.width_s);  % 0 = blanked

% NCO of the DDC: the real line is multiplied by cos(2*pi*fIF*t), the imaginary line by -sin(2*pi*fIF*t)
sys.ddc.cos      = single(cos(2*pi*sys.fIF*sys.tAdc));
sys.ddc.minusSin = single(-sin(2*pi*sys.fIF*sys.tAdc));

% DDC output samples: their time (a valid output is centred on its input window, linear-phase
% FIR) and the first/last ADC sample each of them depends on
t1 = sys.tAdc(1:D1:sys.nAdc-L1+1) + (L1-1)/2/fs;
sys.tDdc = t1(1:D2:numel(t1)-L2+1) + (L2-1)/2/fs1;
sys.nDdc = numel(sys.tDdc);
k = (0:sys.nDdc-1)';
firstAdc = k*D2*D1 + 1;
lastAdc  = (k*D2 + L2 - 1)*D1 + L1;

%% Receive interval, DDC segment, matched filter and delay / range axes of every sub-pulse
rxInterval.short = [sys.rxStart_s, wl.tTx_s];                 % [start end) of the receive interval
rxInterval.long  = [wl.tTx_s + wl.width_s, inf];
for i = 1:2
    wv = sys.wave.(names{i});
    % DDC samples whose whole input window lies inside the receive interval of this sub-pulse
    iv = rxInterval.(names{i});
    wv.ddcIdx = find(sys.tAdc(firstAdc) >= iv(1) & sys.tAdc(lastAdc) < iv(2)).';
    % Matched filter (Taylor window for LFM), unit noise gain
    Lp  = round(wv.width_s * sys.fsOut);                      % matched filter length [samples]
    ref = waveform_envelope(wv, ((0:Lp-1)' + 0.5) / sys.fsOut);  % mid-sample times (not on chip edges)
    if strcmpi(wv.type, 'lfm')
        win = taylor_window(Lp, cfg.proc.taylorNbar, cfg.proc.taylorSll_dB);
    else
        win = ones(Lp, 1);
    end
    h = conj(flipud(ref .* win));
    wv.mfTaps = (h / norm(h)).';                              % row vector
    % Output k of the valid convolution correlates the reference with the
    % segment starting at sample k. The reference holds the pulse at the
    % mid-sample times, so output k peaks for an echo that starts half a
    % sample before tDdc(k).
    wv.nRange  = numel(wv.ddcIdx) - Lp + 1;
    wv.delay_s = sys.tDdc(wv.ddcIdx(1:wv.nRange)).' - 0.5/sys.fsOut - wv.tTx_s;   % delay axis [s]
    wv.range_m = c * wv.delay_s / 2;                              % range axis [m]
    sys.wave.(names{i}) = wv;
end
if sys.wave.short.range_m(end) < sys.wave.long.range_m(1)
    warning('Gap in the range coverage: short pulse up to %.1f km, long pulse from %.1f km.', ...
        sys.wave.short.range_m(end)/1e3, sys.wave.long.range_m(1)/1e3);
end

%% Digital beamforming (disabled with a single receive element)
sys.dbf.enabled = sys.nCh > 1;
if sys.dbf.enabled
    dbfF = load_input(cfg, 'dbf', struct('fc', cfg.fc, 'nElements', sys.nCh, 'tilt_deg', cfg.ant.tilt_deg));
    if size(dbfF.W, 2) ~= sys.nCh
        error('The DBF matrix has %d columns but cfg.ant.nElements = %d.', size(dbfF.W, 2), sys.nCh);
    end
    sys.dbf.W = dbfF.W;
    sys.dbf.beamEl_deg = dbfF.beamEl_deg(:);
    sys.dbf.beamU = sind(sys.dbf.beamEl_deg - sys.ant.tilt_deg);
else
    sys.dbf.W = [];
    sys.dbf.beamEl_deg = NaN;                                 % no elevation measurement
    sys.dbf.beamU = NaN;
end
sys.nBeams = numel(sys.dbf.beamEl_deg);

% All clutter arrives from (about) the same elevation: one steering vector
sys.clutterSteer = steering_vector(sys, 0, cfg.clutter.elevation_deg);

%% CPI buffer and Doppler processing
% CPI = time a target stays inside the 3 dB azimuth beam; the buffer holds
% the pulses that arrive during one CPI.
sys.cpi.time_s = sys.ant.azBeamwidth_deg / sys.ant.rate_dps;
if isempty(cfg.proc.cpiLength)
    sys.cpi.length = round(sys.cpi.time_s / mean(sys.pri));
else
    sys.cpi.length = cfg.proc.cpiLength;
end
nTaps = size(sys.mtiCoeffs, 2);
nOut  = sys.cpi.length - (nTaps - 1);                         % MTI output samples
if sys.staggered
    nRows  = sys.nStagger;                                    % one row per stagger position
    rowLen = floor(nOut / nRows);
    rowPri = sum(sys.pri);                                    % sampling period inside a row
else
    nRows  = 1;                                               % one FFT over the whole CPI
    rowLen = nOut;
    rowPri = sys.pri;
end
nfft = cfg.proc.nfft;
if isempty(nfft) && sys.staggered
    nfft = 2^nextpow2(rowLen);                                % 28 samples per row -> 32
elseif isempty(nfft)
    nfft = 2^nextpow2(sys.cpi.length);                        % 114 pulses -> 128
end
rowLen = min(rowLen, nfft);                                   % a shorter FFT keeps the newest samples
sys.cpi.nRows  = nRows;
sys.cpi.rowLen = rowLen;
w = hamming_window(rowLen);
sys.dop.window = w / norm(w);
sys.dop.nfft   = nfft;
sys.dop.rowPri_s = rowPri;
sys.dop.f_Hz   = ((0:nfft-1) - floor(nfft/2)) / (nfft*rowPri); % after fftshift
sys.dop.v_mps  = -sys.dop.f_Hz * sys.lambda / 2;              % radial velocity (+ = away)

% Noise power of every Doppler bin for white noise at the MTI input. The MTI
% colours the noise (e.g. +4 dB near PRF/2 for [1 -2 1]); doppler_process
% divides by this gain so that the CFAR sees a flat noise floor.
n = sys.cpi.length;
sys.dop.noiseGain = 1;
impulses = reshape(eye(n), 1, n, n);                          % one impulse per "range cell"
st = mod(0:n-1, sys.nStagger) + 1;
[pN, Y] = doppler_process(mti_filter(impulses, sys.mtiCoeffs, st), sys);
sys.dop.noiseGain = sum(pN, 2);
% Effective number of independent looks of the non-coherent integration.
% De-interleaved rows (staggered PRF) are built from overlapping MTI outputs,
% so their noise is correlated. The tail of the averaged noise power decays
% with the largest eigenvalue of the row covariance, so the number of looks
% is trace(C) / max(eig(C)): 1 for a single FFT, 1.5 for the default stagger.
nLooks = zeros(1, nfft);
for k = 1:nfft
    Yk = reshape(Y(1, :, k, :), n, nRows);
    C  = Yk' * Yk;                                            % noise covariance of the rows
    nLooks(k) = real(trace(C)) / max(real(eig(C)));
end
% MTI notch: Doppler bins around zero (clutter residue) that the CFAR does not test.
% The residue is as wide as the Doppler window mainlobe: +-2*nfft/rowLen bins for Hamming.
sys.dop.notchBins = round(cfg.cfar.notchWidth * 2 * nfft / rowLen);
zeroBin = floor(nfft/2) + 1;
tested = abs((1:nfft) - zeroBin) > sys.dop.notchBins;
sys.dop.nLooks = min(nLooks(tested));
% Median / mean of the noise power (used to estimate the noise floor of a beam)
L = sys.dop.nLooks;
sys.dop.medianToMean = fzero(@(x) gammainc(L*x, L) - 0.5, [0.1 2]);

%% CFAR
sys.cfar = cfg.cfar;
nOuter = (2*(cfg.cfar.guardRange + cfg.cfar.trainRange) + 1) * (2*(cfg.cfar.guardDoppler + cfg.cfar.trainDoppler) + 1);
nGuard = (2*cfg.cfar.guardRange + 1) * (2*cfg.cfar.guardDoppler + 1);
sys.cfar.nTrain = nOuter - nGuard;
sys.cfar.alpha  = cfar_alpha(cfg.cfar.pfa, sys.cfar.nTrain, sys.dop.nLooks);
if 2*(cfg.cfar.guardDoppler + cfg.cfar.trainDoppler) + 1 > nfft
    error('CFAR Doppler window is larger than the FFT length.');
end
end

function bw = beamwidth_3dB(az_deg, gain_dB)
% 3 dB width of a pattern, interpolated between the table points.
g = gain_dB - max(gain_dB);
[~, iMax] = max(g);
iL = find(g(1:iMax) < -3.0103, 1, 'last');                    % last point below -3 dB on the left
iR = iMax - 1 + find(g(iMax:end) < -3.0103, 1, 'first');      % first point below -3 dB on the right
azL = interp1(g(iL:iL+1), az_deg(iL:iL+1), -3.0103);
azR = interp1(g(iR-1:iR), az_deg(iR-1:iR), -3.0103);
bw = azR - azL;
end

function data = load_input(cfg, which, madeFor)
% Load input file cfg.files.(which). A missing file, or a default file built by
% make_input_files for other configuration values (its 'madeFor'), is rebuilt.
file = cfg.files.(which);
if exist(file, 'file') ~= 2
    make_input_files(cfg, which);
end
data = load(file);
if isfield(data, 'madeFor') && ~isequal(data.madeFor, madeFor)
    make_input_files(cfg, which);
    data = load(file);
end
end
