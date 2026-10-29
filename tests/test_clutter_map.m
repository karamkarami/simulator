function ok = test_clutter_map(receiverBin)
%TEST_CLUTTER_MAP  Clutter map instead of the MTI filter (cfg.clutterRejection.method = 'map').
%   1. Threshold factors (sys.cmap.alpha): with exponentially distributed
%      powers the map test has the design false alarm probability after 1,
%      2, 4, 8 and 40 averaged scans (Monte Carlo at Pfa = 1e-3).
%   2. Map logic on synthetic power maps (clutter_map_detect.m,
%      clutter_map_fold.m): no detection in an azimuth cell without a map;
%      the band powers of a pass enter the map only when the antenna leaves
%      the cell (a target does not raise its own threshold); the map is the
%      mean of the first scans, then exponential averaging; a band cell is a
%      detection when it passed the CFAR and exceeds alpha(n) times the map;
%      the bins outside the band keep the CFAR decision.
%   3. End to end over 6 scans of a sector with ground clutter (one element):
%      target A flies radially at a blind speed (3 x lambda*PRF/2: its
%      Doppler folds to 0 Hz), target B at another speed. With the MTI only
%      B is detected (A is in the notch). With the clutter map B is detected
%      in every scan and A once the cells it flies into have been mapped in
%      earlier scans; no false plot. When the C receiver is built it gets the
%      DDC stream of the clutter-map run and must give the same plots.

root = fileparts(fileparts(mfilename('fullpath')));
addpath(root, fullfile(root, 'src'));
if nargin < 1 || isempty(receiverBin)
    receiverBin = fullfile(root, 'radar_receiver', 'build', 'radar_receiver');
end
ok = true;
cfg = radar_config();
cfg.sim.verbose    = false;
cfg.display.enable = false;
cfg.io.saveRaw     = false;
cfg.ant.nElements  = 1;                                % (quick: no DBF)
cfg.rx.maxRange_m  = 50e3;
cfg.clutterRejection.method = 'map';

%% 1) threshold factors against a Monte Carlo of exponential powers
c = cfg;
c.clutterRejection.mapPfa = 1e-3;
sys = setup_radar(c);
w = sys.cmap.weight;
rng(5);
N = 2e6;
for n = [1 2 4 8 40]
    m = zeros(N, 1);
    for k = 1:n                                        % the map after n scans (clutter_map_fold.m)
        m = m + max(w, 1/k) * (-log(rand(N, 1)) - m);
    end
    pfa = mean(-log(rand(N, 1)) > sys.cmap.alpha(n) * m);
    ok = report(sprintf('map after %2d scan(s): threshold %5.2f dB, false alarms %.2e (design 1e-3)', n, ...
                        10*log10(sys.cmap.alpha(n)), pfa), abs(pfa/1e-3 - 1) < 0.1) && ok;
end

%% 2) map logic on synthetic power maps
sys = setup_radar(cfg);
wv = sys.wave.long;
band = sys.cmap.bins;
other = setdiff(1:sys.dop.nfft, band);
cm = clutter_map_init(sys, wv);
P = ones(sys.nBeams, wv.nRange, sys.dop.nfft, 'single');
all1 = true(size(P));
az = @(cell, frac) (cell - 1 + frac) * sys.cmap.azCell_deg;     % an azimuth inside map cell 'cell'
[d, cm] = clutter_map_detect(cm, all1, P, az(11, 0.2), sys);    % scan 1, cell 11
ok = report('scan 1: no detection in the band without a map, the other bins keep the CFAR decision', ...
            ~any(any(any(d(:, :, band)))) && all(all(all(d(:, :, other))))) && ok;
P2 = P;  P2(1, 50, band(3)) = 1e9;                     % a strong target in this pass
[~, cm] = clutter_map_detect(cm, all1, P2, az(11, 0.6), sys);
ok = report('the pass over a cell does not enter the map while the antenna is in it', ...
            cm.count(11) == 0 && all(cm.M(:) == 0) && cm.n == 2) && ok;
[~, cm] = clutter_map_detect(cm, all1, P, az(12, 0.1), sys);    % the antenna enters cell 12
mean11 = (1 + 1e9) / 2;                                % mean of the two updates of the pass
M1 = cm.M(:, :, :, 11);
ok = report('leaving the cell: its mean band power is the map (first scan)', cm.count(11) == 1 ...
            && M1(1, 50, 3) == single(mean11) && nnz(M1 ~= 1) == 1) && ok;
cm = clutter_map_fold(cm, sys);                        % (end of the scan)
P3 = P;  P3(1, 20, band(2)) = 1.01 * sys.cmap.alpha(1);          % scan 2: above alpha(1) x map
P3(1, 21, band(2)) = 0.99 * sys.cmap.alpha(1);
isDet3 = all1;  isDet3(1, 22, band(2)) = false;        % the CFAR did not pass this cell
P3(1, 22, band(2)) = 2 * sys.cmap.alpha(1);
[d, cm] = clutter_map_detect(cm, isDet3, P3, az(11, 0.5), sys);
dBand = d(:, :, band);
ok = report(sprintf(['scan 2: band detections only above alpha(1) = %.1f dB x map and with the CFAR: ' ...
                     '%d detections'], 10*log10(sys.cmap.alpha(1)), nnz(dBand)), ...
            nnz(dBand) == 1 && d(1, 20, band(2))) && ok;
cm = clutter_map_fold(cm, sys);
M2 = M1 + single(1/2) * (single(double(P3(:, :, band))) - M1);
ok = report('scan 2 enters the map with weight 1/2 (mean of two scans)', cm.count(11) == 2 ...
            && isequal(cm.M(:, :, :, 11), M2)) && ok;
for s = 3:60                                           % many scans of power 3
    [~, cm] = clutter_map_detect(cm, all1, 3*P, az(11, 0.5), sys);
    cm = clutter_map_fold(cm, sys);
end
ok = report(sprintf('after 60 scans of power 3: map %.4f, scans counted %d (threshold %.1f dB)', ...
                    cm.M(1, 7, 1, 11), cm.count(11), 10*log10(sys.cmap.alpha(end))), ...
            abs(cm.M(1, 7, 1, 11) - 3) < 1e-3 && cm.count(11) == numel(sys.cmap.alpha)) && ok;

%% 3) end to end: a target at a blind speed and one at another speed, MTI and clutter map
vBlind = sys.lambda / (2 * sys.pri(1));                % first blind speed
cfg.sim.nScans   = 6;
cfg.sim.swerling = 0;
cfg.targets = [make_target('A', 30, 40e3, 60, 3000, -3*vBlind, 0), ...    % Doppler folds to 0 Hz
               make_target('B', 10, 35e3, 64, 3000, -120, 0)];
s = setup_radar(cfg);
sectors = target_sectors(s, cfg.targets);
cfg.sim.sectors_deg = [min(sectors(:, 1)), max(sectors(:, 2))];
cfg.io.outFolder = fullfile(root, 'output', 'test_clutter_map');
cMti = cfg;
cMti.clutterRejection.method = 'mti';
resMti = run_radar_sim(cMti);
useC = isunix && exist(receiverBin, 'file');
if useC                                                % the clutter-map run streams to the C receiver
    cfg.udp.enable = true;
    cfg.udp.port   = 50330;
    resMap = c_receiver_live(cfg, receiverBin, cfg.io.outFolder);
else
    resMap = run_radar_sim(cfg);
end
names = {'MTI', 'clutter map'};
res = {resMti, resMap};
for i = 1:2
    [scansA, scansB, nFalse] = plot_scans(res{i}, cfg.targets);
    if i == 1
        good = isempty(scansA) && isequal(scansB, 1:6) && nFalse == 0;
    else
        good = numel(scansA) >= 2 && ~any(scansA == 1) && isequal(scansB, 1:6) && nFalse == 0;
    end
    ok = report(sprintf('%-11s: A (blind speed) plotted in scans [%s], B in scans [%s], %d false plot(s)', ...
                        names{i}, num2str(scansA), num2str(scansB), nFalse), good) && ok;
end
if useC
    f = c_receiver_files(cfg.io.outFolder);
    cmp = c_receiver_compare(c_receiver_reference(resMap), f.plotsLive);
    ok = report(sprintf('C receiver with the clutter map: the %d MATLAB plots', numel(resMap.plots)), cmp.ok) && ok;
end
fprintf('Clutter map: %s\n', ternary(ok, 'PASS', 'FAIL'));
end

function [scansA, scansB, nFalse] = plot_scans(res, targets)
% Scans in which targets A and B have a plot within 300 m and 2 deg; plots elsewhere.
scansA = [];  scansB = [];  nFalse = 0;
scanTime = 60 / res.sys.cfg.ant.rpm;
for p = res.plots
    hit = 0;
    for j = 1:2
        [R, az] = target_geometry(targets(j), p.time_s, res.sys.ant.height_m);
        if abs(p.range_m - R) < 300 && abs(wrap180(p.az_deg - az)) < 2
            hit = j;
        end
    end
    scan = floor(p.time_s / scanTime) + 1;
    if hit == 1, scansA(end+1) = scan; end %#ok<AGROW>
    if hit == 2, scansB(end+1) = scan; end %#ok<AGROW>
    nFalse = nFalse + (hit == 0);
end
scansA = unique(scansA);
scansB = unique(scansB);
end

function ok = report(text, ok)
fprintf('  %-110s %s\n', text, ternary(ok, 'ok', 'FAIL'));
end

function s = ternary(c, a, b)
if c, s = a; else, s = b; end
end
