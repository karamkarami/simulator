function ok = test_snr_budget(prfMode, pulseName, waveType)
%TEST_SNR_BUDGET  Compare the expected SNR of every processing block with a measurement.
%
%   ok = TEST_SNR_BUDGET(prfMode, pulseName) places one target in the middle
%   of the coverage of sub-pulse pulseName ('short' or 'long'), sweeps the
%   antenna over it for one CPI and runs the processing chain twice:
%     - signal only (no thermal noise, no ADC noise, no quantization)
%     - noise only  (no target)
%   The SNR measured after the ADC, DDC, DBF, matched filter, MTI and Doppler
%   blocks is compared with snr_budget.m. prfMode is 'fixed' or 'staggered',
%   waveType is 'lfm' (default) or 'phasecode'.

if nargin < 1, prfMode = 'fixed'; end
if nargin < 2, pulseName = 'long'; end
if nargin < 3, waveType = 'lfm'; end
root = fileparts(fileparts(mfilename('fullpath')));
addpath(root, fullfile(root, 'src'));
rng(7);

cfg = radar_config();
cfg.prf.mode = prfMode;
cfg.wave.(pulseName).type = waveType;
cfg.rx.maxRange_m = 50e3;
cfg.clutter.enable = false;
sys = setup_radar(cfg);
wv  = sys.wave.(pulseName);
c   = cfg.const.c;

%% Target: on a range cell at the CPI centre, broadside of a beam, moving radially
n     = sys.cpi.length;
st    = mod(0:n-1, sys.nStagger) + 1;
pri   = sys.pri(st(1:end-1));
tp    = [0, cumsum(pri(:).')];                       % pulse times
tMid  = mean(tp);
azTgt = 100;                                         % target azimuth [deg]
azAnt = azTgt + sys.ant.rate_dps*(tp - tMid);        % antenna sweeps across the target
R0    = wv.range_m(round(0.6*wv.nRange));            % range at the CPI centre (cell centre)
el    = 6;                                           % target elevation [deg]
vr    = -180;                                        % radial velocity [m/s]
ground = R0*cosd(el);
tgt = make_target('S', 1, ground, azTgt, cfg.ant.height_m + R0*sind(el), 0, 0);
tgt.vel_mps = vr*(tgt.pos0_m - [0 0 cfg.ant.height_m])/R0;   % purely radial motion
tgt.pos0_m  = tgt.pos0_m - tgt.vel_mps*tMid;                 % at R0 at the CPI centre

budget = snr_budget(sys, tgt, tMid);

%% Run one CPI twice: signal only and noise only
cfgS = cfg;  cfgS.rx.addThermalNoise = false;  cfgS.adc.addAdcNoise = false;
sysS = sys;  sysS.cfg = cfgS;
[rpS, xS] = run_cpi(sysS, tgt, tp, azAnt, wv, false);
[rpN, xN] = run_cpi(sys, [], tp, azAnt, wv, true);

%% Measure the SNR after every block (centre pulse of the CPI)
tau = 2*R0/c;
inEcho = @(t) t >= wv.tTx_s + tau + 0.25*wv.width_s & t <= wv.tTx_s + tau + 0.75*wv.width_s;

m.adc = mean(double(xS.adc(inEcho(sys.tAdc), 1)).^2) / var(double(xN.adc(:, 1)));
m.ddc = mean(abs(xS.ddc(inEcho(sys.tDdc(1:size(xS.ddc, 1))), 1)).^2) / mean(abs(xN.ddc(:, 1)).^2);
beam  = budget.beam;
m.dbf = mean(abs(xS.dbf(beam, inEcho(sys.tDdc(1:size(xS.dbf, 2))))).^2) / mean(abs(xN.dbf(beam, :)).^2);
[pk, r] = max(abs(xS.mf(beam, :)).^2);
m.mf  = pk / mean(abs(xN.mf(beam, :)).^2);

yS = mti_filter(rpS, sys.mtiCoeffs, st);
yN = mti_filter(rpN, sys.mtiCoeffs, st);
mid = round(size(yS, 3)/2) + (-2*sys.nStagger:2*sys.nStagger);
noiseMti = abs(yN(beam, :, :)).^2;
m.mti = mean(abs(yS(beam, r, mid)).^2) / mean(noiseMti(:));
pS = doppler_process(yS, sys);
pN = doppler_process(yN, sys);
[pk, kd] = max(pS(beam, r, :));
m.doppler = pk / mean(pN(beam, :, kd));

%% Compare
blocks = {'adc', 'ddc', 'dbf', 'mf', 'mti', 'doppler'};
tol    = [0.5, 0.5, 0.5, 1.0, 1.0, 1.5];               % [dB]
fprintf('\nSNR budget check: %s PRF, %s pulse (%s), R = %.1f km, el = %.1f deg, vr = %.0f m/s, beam %d\n', ...
    prfMode, pulseName, waveType, R0/1e3, el, vr, beam);
fprintf('%-8s %10s %10s %8s\n', 'Block', 'Expected', 'Measured', 'Diff');
ok = true;
for i = 1:numel(blocks)
    e = budget.(blocks{i});
    v = 10*log10(m.(blocks{i}));
    pass = abs(v - e) <= tol(i);
    ok = ok && pass;
    fprintf('%-8s %10.2f %10.2f %8.2f %s\n', blocks{i}, e, v, v - e, ternary(pass, '', '  <-- FAIL'));
end
fprintf('Result: %s\n', ternary(ok, 'PASS', 'FAIL'));
end

function [rpCpi, last] = run_cpi(sys, targets, tp, azAnt, wv, quantize)
% Processing chain up to the matched filter for every pulse of one CPI.
n = numel(tp);
rpCpi = complex(zeros(sys.nBeams, wv.nRange, n, 'single'));
k0 = ceil(n/2);
for k = 1:n
    x = generate_rx_signal(sys, targets, tp(k), azAnt(k), []);
    if quantize
        x = adc_sample(x, sys);
    end
    bb = ddc(x, sys, wv);
    bm = digital_beamforming(bb, sys.dbf.W);
    rp = matched_filter(bm, wv);
    rpCpi(:, :, k) = rp;
    if k == k0                                       % keep the centre pulse of every block
        last.adc = x;  last.ddc = bb;  last.dbf = bm;  last.mf = rp;
    end
end
end

function s = ternary(cond, a, b)
if cond, s = a; else, s = b; end
end
