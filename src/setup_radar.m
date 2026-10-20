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
sys.fLO    = cfg.fc - cfg.adc.fIFcenter_Hz;                   % receiver local oscillator [Hz]

%% Input files (the defaults are created when a file is missing)
files = {cfg.files.antenna, cfg.files.ddc, cfg.files.dbf, cfg.files.stagger};
if ~all(cellfun(@(f) exist(f, 'file') == 2, files))
    make_input_files(cfg);
end
ant  = load(cfg.files.antenna);
ddcF = load(cfg.files.ddc);
dbfF = load(cfg.files.dbf);

%% Antenna
sys.ant = ant;
sys.ant.height_m = cfg.ant.height_m;
sys.ant.tilt_deg = cfg.ant.tilt_deg;
sys.ant.rate_dps = cfg.ant.rpm * 360 / 60;                    % rotation rate [deg/s]
inBeam = ant.az_deg(ant.azGain_dB >= max(ant.azGain_dB) - 3.0103);
sys.ant.azBeamwidth_deg = max(inBeam) - min(inBeam);          % 3 dB azimuth beamwidth
sys.ant.azSpan_deg = 3 * sys.ant.azBeamwidth_deg;             % echoes outside +-azSpan are ignored
sys.nCh = ant.nElements;

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
        st = load(cfg.files.stagger);
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

%% DDC filters
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

%% Timing inside one PRI: long pulse, short pulse, then the receive window
order = {'long', 'short'};
tTx = 0;
for i = 1:numel(order)
    wv = cfg.wave.(order{i});
    wv.name  = order{i};
    wv.tTx_s = tTx;                                           % transmit start in the PRI [s]
    wv.fRF_Hz   = cfg.fc + wv.fIF_Hz - cfg.adc.fIFcenter_Hz;   % carrier of this sub-pulse [Hz]
    wv.lambda_m = c / wv.fRF_Hz;
    tTx = tTx + wv.width_s;
    sys.wave.(order{i}) = wv;
end
sys.rxStart_s = ceil(tTx*fs) / fs;
filterSpan = (L1-1)/fs + (L2-1)/fs1;                          % time lost by the valid convolutions
rxEnd = sys.wave.long.tTx_s + 2*cfg.rx.maxRange_m/c + sys.wave.long.width_s + filterSpan + 2/sys.fsOut;
if rxEnd > min(sys.pri)
    error('Receive window (%.1f us) does not fit in the shortest PRI; reduce cfg.rx.maxRange_m.', rxEnd*1e6);
end
sys.nAdc = round((rxEnd - sys.rxStart_s) * fs);
sys.tAdc = sys.rxStart_s + (0:sys.nAdc-1)' / fs;             % ADC sample times in the PRI [s]

% DDC output times: a valid output sample is centred on its input window (linear-phase FIR)
t1 = sys.tAdc(1:D1:sys.nAdc-L1+1) + (L1-1)/2/fs;
sys.tDdc = t1(1:D2:numel(t1)-L2+1) + (L2-1)/2/fs1;
sys.nDdc = numel(sys.tDdc);

%% Matched filters and the delay / range axis of every sub-pulse
for i = 1:numel(order)
    wv = sys.wave.(order{i});
    Lp   = round(wv.width_s * sys.fsOut);                     % matched filter length [samples]
    ref  = waveform_envelope(wv, (0:Lp-1)' / sys.fsOut);
    if strcmpi(wv.type, 'lfm')
        win = taylor_window(Lp, cfg.proc.taylorNbar, cfg.proc.taylorSll_dB);
    else
        win = ones(Lp, 1);
    end
    h = conj(flipud(ref .* win));
    wv.mfTaps = (h / norm(h)).';                              % row vector, unit noise gain
    % Output k of the valid convolution correlates the reference with the
    % input starting at sample k, i.e. an echo that starts at tDdc(k).
    nMf   = sys.nDdc - Lp + 1;
    delay = sys.tDdc(1:nMf) - wv.tTx_s;
    range = c * delay / 2;
    keep  = find(range >= wv.procRange_m(1) & range <= wv.procRange_m(2));
    wv.keepIdx = keep(:).';
    wv.delay_s = delay(keep).';                               % delay axis [s]
    wv.range_m = range(keep).';                               % range axis [m]
    wv.nRange  = numel(keep);
    % ADC samples needed to compute the last kept range cell (the DDC skips the rest)
    nDdcUsed = keep(end) + Lp - 1;
    wv.nAdcUsed = ((nDdcUsed - 1)*D2 + L2 - 1)*D1 + L1;
    sys.wave.(order{i}) = wv;
end

%% Digital beamforming
sys.dbf.W = dbfF.W;
sys.dbf.beamEl_deg = dbfF.beamEl_deg(:);
sys.dbf.beamU = sind(sys.dbf.beamEl_deg - sys.ant.tilt_deg);
sys.nBeams = size(dbfF.W, 1);
if size(dbfF.W, 2) ~= sys.nCh
    error('DBF matrix has %d columns but the antenna has %d elements.', size(dbfF.W, 2), sys.nCh);
end

% All clutter arrives from (about) the same elevation: one steering vector
sys.clutterSteer = steering_vector(sys, 0, cfg.clutter.elevation_deg);

%% CPI buffer and Doppler processing
nfft = cfg.proc.nfft;
nTaps = size(sys.mtiCoeffs, 2);
if isempty(cfg.proc.cpiLength)
    dwell = sys.ant.azBeamwidth_deg / sys.ant.rate_dps;       % time on target [s]
    nCpi  = floor(dwell / mean(sys.pri));
else
    nCpi = cfg.proc.cpiLength;
end
nOut = nCpi - (nTaps - 1);                                    % MTI output samples
if sys.staggered
    nRows  = sys.nStagger;                                    % one row per stagger position
    rowLen = min(floor(nOut / nRows), nfft);
    rowPri = sum(sys.pri);                                    % sampling period inside a row
else
    nRows  = max(1, floor(nOut / nfft));                      % consecutive blocks of nfft pulses
    rowLen = min(nOut, nfft);
    rowPri = sys.pri;
end
sys.cpi.nRows  = nRows;
sys.cpi.rowLen = rowLen;
sys.cpi.length = nRows*rowLen + nTaps - 1;                    % pulses held by a CPI buffer
w = hamming_window(rowLen);
sys.dop.window = w / norm(w);
sys.dop.nfft   = nfft;
sys.dop.rowPri_s = rowPri;
sys.dop.f_Hz   = ((0:nfft-1) - floor(nfft/2)) / (nfft*rowPri); % after fftshift
sys.dop.v_mps  = -sys.dop.f_Hz * sys.lambda / 2;              % radial velocity at fc (+ = away)
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
% is trace(C) / max(eig(C)): 3 for 3 independent blocks, 1.5 for the default stagger.
nLooks = zeros(1, nfft);
for k = 1:nfft
    Yk = reshape(Y(1, :, k, :), n, nRows);
    C  = Yk' * Yk;                                            % noise covariance of the rows
    nLooks(k) = real(trace(C)) / max(real(eig(C)));
end
zeroBin = floor(nfft/2) + 1;
tested = abs((1:nfft) - zeroBin) > cfg.cfar.notchBins;        % bins outside the MTI notch
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
