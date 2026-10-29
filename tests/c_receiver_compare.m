function cmp = c_receiver_compare(ref, plotsFile, tracksFile)
%C_RECEIVER_COMPARE  Plots and track updates of the C receiver against those of MATLAB.
%   cmp = C_RECEIVER_COMPARE(ref, plotsFile) matches every MATLAB plot
%   (ref from c_receiver_reference.m) with a plot of the C receiver (its
%   --plots-csv file): same sub-pulse, beam/range/Doppler cell +-1 (Doppler
%   circular), azimuth within 0.5 deg, and nearly the same range, azimuth,
%   elevation, radial velocity and SNR. The C receiver must have no other
%   plot. cmp.ok is the result, cmp.maxDiff the largest differences.
%   cmp = C_RECEIVER_COMPARE(ref, plotsFile, tracksFile) also compares the
%   track updates (--tracks-csv) with the MATLAB track log: the same track
%   numbers, update by update, with positions within 20 m (cmp.tracksOk).
cmp = compare_plots(ref.plots, read_plots_csv(plotsFile), ref.axes);
if nargin > 2 && ~isempty(tracksFile)
    cmp.tracksOk = compare_tracks(ref.trackLog, tracksFile);
    cmp.ok = cmp.ok && cmp.tracksOk;
end
end

function cmp = compare_plots(pm, pc, ax)
% Matches every MATLAB plot with a C plot (same sub-pulse, beam/range/Doppler
% cell +-1 (Doppler circular), azimuth within 0.5 deg) and checks the values.
tol.range_m = 1.01 * ax.rangeCell_m;
tol.az_deg = 0.2;  tol.el_deg = 0.2;  tol.snr_dB = 0.5;
tol.vr_mps = 1.01 * ax.velocityBin_mps;
nfft = ax.nfft;
usedC = false(1, numel(pc));
cmp.ok = true;
cmp.maxDiff = struct('range_m', 0, 'az_deg', 0, 'el_deg', 0, 'vr_mps', 0, 'snr_dB', 0, 'hits', 0);
fprintf('\n%5s %6s | %9s %9s | %7s %7s | %6s %6s | %7s %7s | %5s %5s | %4s %4s | %s\n', 'pulse', 't[s]', ...
    'R MATLAB', 'R C', 'az M', 'az C', 'el M', 'el C', 'vr M', 'vr C', 'SNR M', 'SNR C', 'hitM', 'hitC', 'result');
for i = 1:numel(pm)
    m = pm(i);
    best = 0;
    for j = find(~usedC)
        c = pc(j);
        dD = abs(c.dIdx + 1 - m.dIdx);  dD = min(dD, nfft - dD);
        if strcmp(c.pulse, m.pulse) && abs(c.beam + 1 - m.beam) <= 1 && abs(c.rIdx + 1 - m.rIdx) <= 1 ...
                && dD <= 1 && abs(wrap180(c.az_deg - m.az_deg)) < 0.5
            if best == 0 || abs(c.time_s - m.time_s) < abs(pc(best).time_s - m.time_s)
                best = j;
            end
        end
    end
    if best == 0
        cmp.ok = false;
        fprintf('%5s %6.2f | %9.1f %9s | no C plot  <-- FAIL\n', m.pulse, m.time_s, m.range_m, '-');
        continue
    end
    usedC(best) = true;
    c = pc(best);
    d.range_m = abs(c.range_m - m.range_m);
    d.az_deg = abs(wrap180(c.az_deg - m.az_deg));
    d.el_deg = abs(c.el_deg - m.el_deg);
    if isnan(c.el_deg) && isnan(m.el_deg), d.el_deg = 0; end
    d.vr_mps = abs(c.vr_mps - m.vr_mps);
    d.snr_dB = abs(c.snr_dB - m.snr_dB);
    d.hits = abs(c.hits - m.hits);
    good = d.range_m <= tol.range_m && d.az_deg <= tol.az_deg && d.el_deg <= tol.el_deg ...
        && d.vr_mps <= tol.vr_mps && d.snr_dB <= tol.snr_dB;
    cmp.ok = cmp.ok && good;
    for f = fieldnames(d).'
        cmp.maxDiff.(f{1}) = max(cmp.maxDiff.(f{1}), d.(f{1}));
    end
    fprintf('%5s %6.2f | %9.1f %9.1f | %7.2f %7.2f | %6.2f %6.2f | %7.2f %7.2f | %5.1f %5.1f | %4d %4d | %s\n', ...
        m.pulse, m.time_s, m.range_m, c.range_m, m.az_deg, c.az_deg, m.el_deg, c.el_deg, m.vr_mps, c.vr_mps, ...
        m.snr_dB, c.snr_dB, m.hits, c.hits, ternary(good, 'ok', '<-- FAIL'));
end
for j = find(~usedC)
    c = pc(j);
    cmp.ok = false;
    fprintf('extra C plot: %s t %.2f s, az %.2f deg, R %.1f m, vr %.2f m/s, SNR %.1f dB  <-- FAIL\n', ...
        c.pulse, c.time_s, c.az_deg, c.range_m, c.vr_mps, c.snr_dB);
end
cmp.nMatlab = numel(pm);
cmp.nC = numel(pc);
fprintf('Plots: MATLAB %d, C %d; largest differences: range %.1f m, az %.3f deg, el %.3f deg, vr %.3f m/s, SNR %.2f dB, hits %d\n', ...
    cmp.nMatlab, cmp.nC, cmp.maxDiff.range_m, cmp.maxDiff.az_deg, cmp.maxDiff.el_deg, cmp.maxDiff.vr_mps, ...
    cmp.maxDiff.snr_dB, cmp.maxDiff.hits);
end

function ok = compare_tracks(logM, fileC)
% MATLAB track log (track_update.m) against the track updates of the C receiver (--tracks-csv).
txt = strtrim(fileread(fileC));
lines = strsplit(txt, {'\n', '\r'});
C = zeros(0, 12);
for i = 2:numel(lines)
    if isempty(strtrim(lines{i})), continue; end
    v = str2double(strsplit(lines{i}, ','));
    C(end+1, :) = v(1:12); %#ok<AGROW>                  % (then the acceleration)
end
C = C(C(:, 3) == 0, :);                                % updates only (not the deletions)
n = min(size(C, 1), size(logM, 1));
sameId = all(C(1:n, 1) == logM(1:n, 2)) && all(C(1:n, 2) == logM(1:n, 3)) && all(C(1:n, 11) == logM(1:n, 10));
dPos = 0;
if n > 0, dPos = max(sqrt(sum((C(1:n, 5:7) - logM(1:n, 4:6)).^2, 2))); end
ok = size(C, 1) == size(logM, 1) && sameId && dPos < 20;
fprintf('Track updates: MATLAB %d, C %d (%d confirmed tracks); same tracks: %s, largest position difference %.2f m: %s\n', ...
    size(logM, 1), size(C, 1), numel(unique(logM(logM(:, 3) == 1, 2))), ternary(sameId, 'yes', 'NO'), dPos, ...
    ternary(ok, 'ok', 'FAIL'));
end

function p = read_plots_csv(fileName)
% Plot list written by the C receiver (--plots-csv), 0-based cell indices.
p = struct('time_s', {}, 'az_deg', {}, 'range_m', {}, 'el_deg', {}, 'vr_mps', {}, 'snr_dB', {}, ...
           'pulse', {}, 'hits', {}, 'beam', {}, 'rIdx', {}, 'dIdx', {}, 'pulseNumber', {});
if ~exist(fileName, 'file'), return; end
lines = strsplit(strtrim(fileread(fileName)), {'\n', '\r'});
for i = 2:numel(lines)                                   % line 1: column names
    if isempty(strtrim(lines{i})), continue; end
    f = strsplit(lines{i}, ',');
    v = str2double(f);
    p(end+1) = struct('time_s', v(1), 'az_deg', v(2), 'range_m', v(3), 'el_deg', v(4), 'vr_mps', v(5), ...
        'snr_dB', v(6), 'pulse', f{7}, 'hits', v(8), 'beam', v(9), 'rIdx', v(10), 'dIdx', v(11), ...
        'pulseNumber', v(12)); %#ok<AGROW>
end
end

function s = ternary(cond, a, b)
if cond, s = a; else, s = b; end
end
