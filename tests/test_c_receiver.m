function [ok, cmp] = test_c_receiver(prfMode, receiverBin, port, targetIdx, replayOnly, cfarMethod)
%TEST_C_RECEIVER  The C receiver (radar_receiver/) against the MATLAB processing.
%
%   [ok, cmp] = TEST_C_RECEIVER(prfMode, receiverBin, port) runs the
%   six-target scenario of test_six_targets ('fixed' or 'staggered' PRF) and
%   streams the DDC output of every PRI over UDP (int16, 1400-byte payloads)
%   to the C receiver, which processes it live and records the datagrams.
%   Every plot of the MATLAB processing must have a C plot in the same
%   beam/range/Doppler cell (+-1) with nearly the same range, azimuth,
%   elevation, radial velocity and SNR, and the C receiver must have no
%   extra plots. The recording is then replayed through the receiver: the
%   plots must be identical to the live ones, and the processing time per
%   PRI is reported.
%
%   receiverBin: the radar_receiver program (default
%   radar_receiver/build/radar_receiver, see radar_receiver/README.md).
%   port: UDP port (default 50000 fixed, 50001 staggered).
%   targetIdx: simulate only the sectors of these targets (default: all six).
%   replayOnly = true: skip the simulation and the live run, and only replay
%   the recording of an earlier run against its saved MATLAB plots (to check
%   a changed receiver quickly).
%   cfarMethod: '2d' (default), 'range' or 'doppler' (cfg.cfar.method).
%   The tracks (EKF, cfg.track) of both programs are compared as well: the
%   same track numbers, update by update, with positions within 20 m.
%   Needs a POSIX shell (Linux, macOS) to start the receiver in the background.

if nargin < 1, prfMode = 'fixed'; end
root = fileparts(fileparts(mfilename('fullpath')));
addpath(root, fullfile(root, 'src'));
if nargin < 2 || isempty(receiverBin)
    receiverBin = fullfile(root, 'radar_receiver', 'build', 'radar_receiver');
end
if nargin < 3 || isempty(port), port = 50000 + strcmp(prfMode, 'staggered'); end
if nargin < 4, targetIdx = []; end
if nargin < 5, replayOnly = false; end
if nargin < 6 || isempty(cfarMethod), cfarMethod = '2d'; end
if ~exist(receiverBin, 'file')
    error('C receiver %s not found: build it first (radar_receiver/README.md)', receiverBin);
end

%% Configuration: the six-target scenario, DDC output streamed over UDP
cfg = radar_config();
cfg.prf.mode        = prfMode;
cfg.cfar.method     = cfarMethod;
cfg.clutter.enable  = true;
cfg.rx.maxRange_m   = 50e3;
cfg.targets         = six_targets_scenario();
cfg.io.outFolder    = fullfile(root, 'output', ['c_receiver_' prfMode '_' cfarMethod]);
cfg.io.saveRaw      = false;
cfg.display.enable  = false;
cfg.udp.enable      = true;
cfg.udp.port        = port;
if isempty(targetIdx), targetIdx = 1:numel(cfg.targets); end
cfg.sim.sectors_deg = target_sectors(setup_radar(cfg), cfg.targets(targetIdx));
out = cfg.io.outFolder;
f = c_receiver_files(out);
if replayOnly
    if ~exist(f.rec, 'file') || ~exist(f.ref, 'file')
        error('No recording of an earlier run in %s', out);
    end
    ref = load(f.ref);
    export_receiver_config(cfg, f.cfg);             % the current format for the current receiver
    ok = true;
    liveDone = false;
else
    if exist(f.ref, 'file'), delete(f.ref); end

    %% Live: start the receiver, run the simulation, wait for the receiver
    results = c_receiver_live(cfg, receiverBin, out);
    ref = c_receiver_reference(results);            % MATLAB plots and the axes needed to compare
    save(f.ref, '-struct', 'ref');

    %% Compare the plots and the tracks
    tracksFile = '';
    if results.sys.track.enable, tracksFile = f.tracksLive; end
    cmp = c_receiver_compare(ref, f.plotsLive, tracksFile);
    ok = cmp.ok;
    liveDone = true;
end

%% Replay of the recording: same plots, processing time
logText = c_receiver_replay(receiverBin, f.cfg, f.rec, f.plotsReplay, f.replayLog);
if liveDone
    same = exist(f.plotsReplay, 'file') && strcmp(fileread(f.plotsReplay), fileread(f.plotsLive));
    fprintf('Replay plots identical to the live plots: %s\n', ternary(same, 'yes', 'NO'));
    ok = ok && same;
else
    cmp = c_receiver_compare(ref, f.plotsReplay);
    ok = cmp.ok;
end
tok = regexp(logText, 'Processing time per PRI: ([\d.]+) ms', 'tokens', 'once');
if ~isempty(tok), cmp.msPerPri = str2double(tok{1}); end
fprintf('Result (%s PRF, C receiver vs MATLAB): %s\n', prfMode, ternary(ok, 'PASS', 'FAIL'));
end

function s = ternary(cond, a, b)
if cond, s = a; else, s = b; end
end
