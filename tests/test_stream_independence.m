function ok = test_stream_independence(receiverBin)
%TEST_STREAM_INDEPENDENCE  The DDC stream depends on the radar only, not on the receiver's processing.
%   1. The simulator (the sender) gives the same ADC record, DDC output and
%      PRI header, bit for bit, whatever the processing settings: CFAR method
%      and window, pulse compression filter (matched, or mismatch filters
%      three times the code length), clutter rejection (MTI or clutter map).
%   2. Every compression filter gives the range cells of the matched filter
%      (samples outside the receive interval count as zeros): with mismatch
%      filters three times the code length a target at the first and at the
%      last range cell of both receive intervals peaks on its cell, also with
%      the default 80 km window (before, the window grew by the filter and
%      did not fit into the PRI).
%   3. When the C receiver is built: the DDC stream of a sender with the 2-D
%      CFAR, matched filters and the MTI (two targets, one element) is
%      recorded and replayed into receivers with other settings. Every PRI is
%      accepted and the plots equal those of the MATLAB processing with the
%      receiver's settings: 1-D Doppler CFAR; 1-D range CFAR with the mismatch
%      filters; clutter map. A receiver with a smaller instrumented range
%      (40 km of the sender's 50 km) accepts the longer record and finds the
%      same plots; one with a larger instrumented range (60 km) needs DDC
%      samples that the stream does not have and ignores every PRI.

root = fileparts(fileparts(mfilename('fullpath')));
addpath(root, fullfile(root, 'src'));
if nargin < 1 || isempty(receiverBin)
    receiverBin = fullfile(root, 'radar_receiver', 'build', 'radar_receiver');
end
ok = true;
out = fullfile(root, 'output', 'test_stream_independence');
if ~exist(out, 'dir'), mkdir(out); end

base = radar_config();
base.sim.verbose    = false;
base.display.enable = false;
base.io.saveRaw     = false;
base.io.outFolder   = out;
base.ant.nElements  = 1;                               % (quick: no DBF)
base.rx.maxRange_m  = 50e3;
base.sim.swerling   = 0;
base.files.phaseCode = fullfile(out, 'phase_code_3x.mat');
for p = {'short', 'long'}
    base.wave.(p{1}).type = 'phasecode';
end
make_codes(base);                                      % P4 codes, mismatch filters 3x the code length
tg = six_targets_scenario();
base.targets = tg([1 3]);                              % 12 km (short pulse) and 35 km (long pulse)

receivers = {'1-D Doppler CFAR',                      @(c) set_cfg(c, 'cfar.method', 'doppler');
             '1-D range CFAR, mismatch filters (3x)', @(c) mismatch(set_cfg(c, 'cfar.method', 'range'));
             'clutter map instead of the MTI',        @(c) set_cfg(c, 'clutterRejection.method', 'map')};

%% 1) the same DDC output for every processing setting
ref = one_pri(base);
for k = 1:size(receivers, 1)
    r = one_pri(receivers{k, 2}(base));
    same = isequal(r.adc, ref.adc) && isequal(r.I, ref.I) && isequal(r.Q, ref.Q) && isequal(r.hdr, ref.hdr);
    ok = report(sprintf('sender with %s: ADC record, DDC output (%d samples) and header identical', ...
                        receivers{k, 1}, size(r.I, 1)), same) && ok;
end

%% 2) mismatch filters three times the code length: the range cells of the matched filter
for maxR = [50e3 80e3]
    c = set_cfg(base, 'rx.maxRange_m', maxR);
    sysM = setup_radar(c);
    sysX = setup_radar(mismatch(c));
    for p = {'short', 'long'}
        wm = sysM.wave.(p{1});  wx = sysX.wave.(p{1});
        sameCells = isequal(wm.range_m, wx.range_m) && sysM.nDdc == sysX.nDdc;
        err = edge_targets(mismatch(c), p{1});
        ok = report(sprintf(['%.0f km, %-5s pulse, %d-tap mismatch filter: the %d cells of the matched filter ' ...
                             '(%.2f ... %.2f km); targets at the first/last cell peak %.2f / %.2f cells off'], ...
                            maxR/1e3, p{1}, numel(wx.mfTaps), wx.nRange, wx.range_m([1 end])/1e3, err), ...
                    sameCells && all(err <= 0.5)) && ok;
    end
end

%% 3) the C receiver with other settings than the sender
if ~(isunix && exist(receiverBin, 'file'))
    fprintf('Stream independence: C receiver not built, its part is skipped\n');
else
    sender = set_cfg(set_cfg(base, 'udp.enable', true), 'udp.port', 50320);
    sender.sim.sectors_deg = target_sectors(setup_radar(sender), sender.targets);
    fS = c_receiver_files(fullfile(out, 'sender'));
    resS = c_receiver_live(sender, receiverBin, fullfile(out, 'sender'));
    cmp = c_receiver_compare(c_receiver_reference(resS), fS.plotsLive);
    ok = report(sprintf('sender settings: C receiver (live) gives the %d MATLAB plots', numel(resS.plots)), ...
                cmp.ok && numel(resS.plots) == 2) && ok;
    nPri = sum(resS.nProcessed);
    for k = 1:size(receivers, 1)
        c = receivers{k, 2}(sender);
        c.udp.enable = false;
        res = run_radar_sim(c);                        % MATLAB: the same DDC data, the receiver's processing
        [cmp, ignored] = replay(receiverBin, c, fS.rec, fullfile(out, sprintf('receiver_%d', k)), res);
        ok = report(sprintf('receiver with %s: %d PRIs ignored, the %d MATLAB plots', receivers{k, 1}, ...
                            ignored, numel(res.plots)), cmp.ok && ignored == 0 && numel(res.plots) >= 2) && ok;
    end
    [cmp, ignored] = replay(receiverBin, set_cfg(sender, 'rx.maxRange_m', 40e3), fS.rec, ...
                            fullfile(out, 'receiver_40km'), resS);
    ok = report(sprintf('receiver with a 40 km window: %d PRIs ignored, the plots of the sender', ignored), ...
                cmp.ok && ignored == 0) && ok;
    [~, ignored, logText] = replay(receiverBin, set_cfg(sender, 'rx.maxRange_m', 60e3), fS.rec, ...
                                   fullfile(out, 'receiver_60km'), []);
    ok = report(sprintf('receiver with a 60 km window: %d of %d PRIs ignored (record too short), no plot', ...
                        ignored, nPri), ignored == nPri && ~isempty(strfind(logText, 'Plots: 0'))) && ok;
end
fprintf('Stream independence: %s\n', ternary(ok, 'PASS', 'FAIL'));
end

function r = one_pri(cfg)
% One PRI with the target, clutter and noise: ADC record, DDC output, header.
cfg.sim.sectors_deg = [];
sys = setup_radar(cfg);
sched = pulse_schedule(sys);
rng(5);
clut = clutter_init(sys);
targets = draw_target_rcs(cfg.targets, 1, 0);
k = 830;                                               % the beam is on target 1
t = sched.time_s(k);  az = sched.az_deg(k);
r.adc = adc_sample(generate_rx_signal(sys, targets, t, az, clutter_signal(clut, sys, az, t)), sys);
[r.I, r.Q] = ddc(r.adc, sys);
r.hdr = pulse_header(sys, sched, k);
end

function err = edge_targets(cfg, pulse)
% Targets without noise at the first and the last range cell of a sub-pulse:
% offset of the compression filter peak from the target cell [cells].
cfg.rx.addThermalNoise = false;  cfg.rx.addLeakage = false;  cfg.adc.addAdcNoise = false;
cfg.adc.quantize = false;        cfg.ddc.fixedPoint = false;  cfg.clutter.enable = false;
sys = setup_radar(cfg);
wv = sys.wave.(pulse);
cells = [1 wv.nRange];
err = zeros(1, 2);
for i = 1:2
    R = wv.range_m(cells(i));
    tgt = draw_target_rcs(make_target('T', 10, sqrt(R^2 - cfg.ant.height_m^2), 0, 0, 0, 0), 1, 0);
    [Rt, ~, ~] = target_geometry(tgt, 0, cfg.ant.height_m);
    [I, Q] = ddc(adc_sample(generate_rx_signal(sys, tgt, 0, 0, []), sys), sys);
    bb = complex(I, Q);
    [~, k] = max(abs(matched_filter(bb(wv.ddcIdx).', wv)));
    err(i) = abs(wv.range_m(k) - Rt) / abs(diff(wv.range_m(1:2)));
end
end

function [cmp, ignored, logText] = replay(receiverBin, cfg, recFile, out, res)
% Replays recFile into the C receiver with the configuration of cfg; compares
% its plots with the MATLAB results res (if given). ignored: PRIs it rejected.
if ~exist(out, 'dir'), mkdir(out); end
f = c_receiver_files(out);
evalc('export_receiver_config(cfg, f.cfg)');
logText = c_receiver_replay(receiverBin, f.cfg, recFile, f.plotsReplay, f.replayLog);
tok = regexp(logText, 'PRIs ignored: (\d+)', 'tokens', 'once');
ignored = 0;
if ~isempty(tok), ignored = str2double(tok{1}); end
cmp.ok = true;
if ~isempty(res)
    cmp = c_receiver_compare(c_receiver_reference(res), f.plotsReplay);
end
end

function make_codes(cfg)
% P4 codes of both sub-pulses with mismatch filters three times the code length.
for p = {'short', 'long'}
    N = round(cfg.wave.(p{1}).width_s / cfg.wave.(p{1}).chip_s);
    m = (0:N-1)';
    code = exp(1j*(pi*m.^2/N - pi*m));                 % as make_input_files.m
    s.([p{1} 'Code']) = code;
    s.([p{1} 'Mismatch']) = mismatch_filter(code, N, 0, 40);   % N chips on each side
end
save(cfg.files.phaseCode, '-struct', 's', '-v7');
end

function c = mismatch(c)
for p = {'short', 'long'}
    c.wave.(p{1}).compression = 'mismatch';
end
end

function c = set_cfg(c, name, value)
parts = strsplit(name, '.');
c = setfield(c, parts{:}, value);
end

function ok = report(text, ok)
fprintf('  %-118s %s\n', text, ternary(ok, 'ok', 'FAIL'));
end

function s = ternary(c, a, b)
if c, s = a; else, s = b; end
end
