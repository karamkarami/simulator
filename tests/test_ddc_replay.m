function ok = test_ddc_replay()
%TEST_DDC_REPLAY  Record the DDC output to files and replay it over UDP with the PRI timing.
%
%   A short run (1 element, one target, no clutter, a 0.5 deg sector, about
%   14 PRIs) saves the DDC output of every PRI (cfg.io.saveDdc) and no raw
%   data (cfg.io.saveRaw = false). The files must hold the same I/Q as the
%   DDC output, replay_ddc_udp.m must send them to localhost bit-exact with
%   their header, and the PRIs must leave on time: median delay below 10 % of
%   the PRI (single late PRIs come from the operating system and are allowed)
%   and, with 'catchUp', the replay must last as long as the recording.

root = fileparts(fileparts(mfilename('fullpath')));
addpath(root, fullfile(root, 'src'));

cfg = radar_config();
cfg.ant.nElements   = 1;
cfg.clutter.enable  = false;
cfg.display.enable  = false;
cfg.sim.verbose     = false;
cfg.sim.sectors_deg = [10 10.5];
cfg.targets         = make_target('A', 10, 20e3, 10.2, 1000, 0, 0);
cfg.io.outFolder    = fullfile(root, 'output', 'test_ddc_replay');
cfg.io.saveRaw      = false;
cfg.io.rawFolder    = fullfile(cfg.io.outFolder, 'raw_pulses');   % must not be created
cfg.io.saveDdc      = true;
cfg.io.ddcFolder    = fullfile(cfg.io.outFolder, 'ddc_pulses');
if exist(cfg.io.rawFolder, 'dir')
    rmdir(cfg.io.rawFolder, 's');
end
if exist(cfg.io.ddcFolder, 'dir')
    delete(fullfile(cfg.io.ddcFolder, 'ddc_*.mat'));
end
res = run_radar_sim(cfg);
files = dir(fullfile(cfg.io.ddcFolder, 'ddc_*.mat'));
filesOk = numel(files) == res.nProcessed && ~exist(cfg.io.rawFolder, 'dir');
d = load(fullfile(cfg.io.ddcFolder, files(1).name));
filesOk = filesOk && isa(d.I, 'int16') && isequal(size(d.I), [res.sys.nDdc, 1]) && isequal(size(d.Q), size(d.I));

port = 50124;
try
    rx = udp_receiver_open(port, 2);
catch err
    fprintf('DDC replay test: SKIPPED (no UDP support: %s)\n', err.message);
    ok = filesOk;
    return
end
% The operating system can stall the sender for milliseconds now and then:
% the timing is checked on the best of up to 3 replays, the data on every one.
dataOk = true;
timingOk = false;
pri = mean(res.sys.pri);
expected = (numel(files) - 1) * pri;
for attempt = 1:3
    stats = replay_ddc_udp(cfg.io.ddcFolder, 'host', '127.0.0.1', 'port', port, 'catchUp', true, 'verbose', false);
    for k = 1:numel(files)
        [hdr, I, Q] = udp_receive_pulse(rx);
        d = load(fullfile(cfg.io.ddcFolder, files(k).name));
        dataOk = dataOk && isequal(I, d.I) && isequal(Q, d.Q) && hdr.pulseNumber == d.hdr.pulseNumber ...
                 && hdr.time_s == d.hdr.time_s && isequal(hdr.shortDdcIdx, d.hdr.shortDdcIdx);
    end
    timingOk = median(stats.late_s) < 0.1*pri && abs(stats.sent_s(end) - expected) < 1e-4 + stats.late_s(end);
    if timingOk, break; end
end
clear rx
ok = filesOk && dataOk && timingOk;
fprintf(['DDC recording and replay (%d PRIs in %.3f ms for %.3f ms of PRIs, delay median %.1f us, max %.1f us, ', ...
         'replay %d): files %s, data %s, timing %s: %s\n'], numel(files), 1e3*stats.sent_s(end), 1e3*expected, ...
        1e6*median(stats.late_s), 1e6*stats.maxLate_s, attempt, ternary(filesOk, 'ok', 'FAIL'), ternary(dataOk, 'ok', 'FAIL'), ...
        ternary(timingOk, 'ok', 'FAIL'), ternary(ok, 'PASS', 'FAIL'));
end

function s = ternary(cond, a, b)
if cond, s = a; else, s = b; end
end
