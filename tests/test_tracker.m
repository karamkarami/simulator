function ok = test_tracker()
%TEST_TRACKER  EKF tracker (track_update.m) on synthetic plots, 3-D and 2-D.
%   The six targets of six_targets_scenario.m are seen once per scan (when
%   the antenna passes them) for 5 scans, with the plot errors of
%   cfg.track (range, azimuth, elevation), plus 8 random false plots per
%   scan. Every target must end with one confirmed track close to its true
%   position and velocity, and no false plot may give a confirmed track.
%   3-D: range, azimuth, elevation; 2-D: range and azimuth (the elevation
%   is unknown, as with a single element, so the track lies at slant range).
%   When the C receiver is built, radar_receiver/build/test_tracker runs the
%   C tracker on the same plots: every track update must match MATLAB.

root = fileparts(fileparts(mfilename('fullpath')));
addpath(root, fullfile(root, 'src'));
ok = true;
for dim = [3 2]
    ok = run_case(dim) && ok;
end
fprintf('EKF tracker: %s\n', ternary(ok, 'PASS', 'FAIL'));
end

function ok = run_case(dim)
cfg = radar_config();
cfg.sim.verbose = false;
cfg.track.dim = sprintf('%dd', dim);
if dim == 2, cfg.ant.nElements = 1; end              % no elevation
sys = setup_radar(cfg);
tgts = six_targets_scenario();
p = sys.track;
rng(7);
tr = track_init(p);
T = p.scanTime_s;
allPlots = zeros(0, 4);                               % plots in processing order (for the C tracker)
rate = sys.ant.rate_dps;
nScans = 5;
for scan = 0:nScans-1
    pl = [];
    for i = 1:numel(tgts)                             % time when the antenna points at the target
        t = scan*T;
        for it = 1:3
            [~, az] = target_geometry(tgts(i), t, cfg.ant.height_m);
            t = scan*T + mod(az - cfg.ant.az0_deg, 360) / rate;
        end
        [R, az, el] = target_geometry(tgts(i), t, cfg.ant.height_m);
        q = make_plot(t, R + cfg.track.sigmaRange_m*randn, az + cfg.track.sigmaAz_deg*randn, ...
                      el + cfg.track.sigmaEl_deg*randn, dim);
        pl = [pl, q]; %#ok<AGROW>
    end
    for f = 1:8                                       % false plots
        az = 360*rand;
        pl = [pl, make_plot(scan*T + az/rate, 5e3 + 45e3*rand, az, 20*rand, dim)]; %#ok<AGROW>
    end
    [~, order] = sort([pl.time_s]);
    for k = order                                     % plots in time order, maintenance every plot
        tr = track_update(tr, pl(k), pl(k).time_s);
        allPlots(end+1, :) = [pl(k).time_s, pl(k).az_deg, pl(k).range_m, pl(k).el_deg]; %#ok<AGROW>
    end
end
tEnd = nScans*T;
tr = track_update(tr, [], tEnd);
rep = track_report(tr);
conf = rep([rep.confirmed]);
ok = numel(conf) == numel(tgts);
fprintf('  %d-D: %d tracks, %d confirmed (expected %d), %d updates\n', dim, numel(rep), numel(conf), ...
        numel(tgts), size(tr.log, 1));
for i = 1:numel(tgts)
    best = inf;  bi = 0;
    for j = 1:numel(conf)
        pos0 = tgts(i).pos0_m + tgts(i).vel_mps*conf(j).t_s - [0 0 cfg.ant.height_m];
        dist = norm(conf(j).pos_m(1:2) - pos0(1:2));
        if dist < best, best = dist; bi = j; end
    end
    if bi == 0, ok = false; continue; end
    tk = conf(bi);
    pos0 = tgts(i).pos0_m + tgts(i).vel_mps*tk.t_s;
    if dim == 2                       % without elevation the track is at slant range (as the plots)
        [R, az] = target_geometry(tgts(i), tk.t_s, cfg.ant.height_m);
        best = norm(tk.pos_m(1:2) - R*[sind(az) cosd(az)]);
    end
    dv = norm(tk.vel_mps(1:2) - tgts(i).vel_mps(1:2));
    dh = NaN;
    if dim == 3, dh = tk.alt_m - pos0(3); end
    good = best < 250 && dv < 25 && (dim == 2 || abs(dh) < 400);
    fprintf('    %s: track %d, position error %5.0f m, velocity error %5.1f m/s, height error %5.0f m %s\n', ...
            tgts(i).name, tk.id, best, dv, dh, ternary(good, 'ok', 'FAIL'));
    ok = ok && good;
end
ok = compare_with_c(allPlots, tr.log, dim, p.gate) && ok;
end

function ok = compare_with_c(plots, log, dim, gate)
% The C tracker on the same plots (when built): it must give the same track updates.
ok = true;
root = fileparts(fileparts(mfilename('fullpath')));
bin = fullfile(root, 'radar_receiver', 'build', 'test_tracker');
if ispc, bin = [bin '.exe']; end
if ~exist(bin, 'file')
    fprintf('    C tracker not built: comparison skipped\n');
    return
end
folder = tempname();
mkdir(folder);
fp = fullfile(folder, 'plots.csv');  fl = fullfile(folder, 'log.csv');
write_csv(fp, plots, '%.17g,%.17g,%.17g,%.17g\n');
write_csv(fl, log, [repmat('%.17g,', 1, 10) '%.17g\n']);
cfgFile = fullfile(root, 'radar_receiver', 'config', 'radar_config.bin');
[status, out] = system(sprintf('"%s" "%s" "%s" "%s" %d %.17g', bin, cfgFile, fp, fl, dim, gate));
fprintf('    C tracker: %s', regexprep(out, '\n(?=.)', '\n    '));
ok = status == 0;
delete(fp);  delete(fl);  rmdir(folder);
end

function write_csv(file, M, fmt)
fid = fopen(file, 'w');
fprintf(fid, fmt, M.');                               % NaN is written as NaN (read by strtod)
fclose(fid);
end

function q = make_plot(t, R, az, el, dim)
if dim == 2, el = NaN; end
q = struct('time_s', t, 'az_deg', mod(az, 360), 'range_m', R, 'el_deg', el, 'vr_mps', 0, 'snr_dB', 20, ...
           'power', 1, 'beam', 1, 'rIdx', 1, 'dIdx', 1, 'pulse', 'long', 'hits', 60);
end

function s = ternary(c, a, b)
if c, s = a; else, s = b; end
end
