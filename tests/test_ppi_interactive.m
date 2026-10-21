function ok = test_ppi_interactive()
%TEST_PPI_INTERACTIVE  Interaction with the PPI: data tips, zoom, pan, filters.
%   Creates the PPI off-screen with a few synthetic plots and drives it
%   through ppi_callbacks.m like the mouse, keyboard and panel would:
%   hover data tips (plot and true target), wheel zoom around the cursor,
%   pan, keys, show/hide of short- and long-pulse plots, SNR and age
%   filters, clear, pinned data tips, labels, reset and PNG export.

root = fileparts(fileparts(mfilename('fullpath')));
addpath(root, fullfile(root, 'src'));
cfg = radar_config();
cfg.display.visible = false;
cfg.targets = make_target('T1', 1, 20e3, 60, 1000, 100, 0);
sys = setup_radar(cfg);
ppi = ppi_init(sys);
fig = ppi.fig;

% Three plots: short at 10 km / 30 deg, long at 40 km / 200 deg (low SNR, old)
p = struct('time_s', {2.0, 0.5, 2.1}, 'az_deg', {30, 200, 95}, 'range_m', {10e3, 40e3, 30e3}, ...
           'el_deg', {NaN, 3, 5}, 'vr_mps', {12.5, -7.25, 3}, 'snr_dB', {25, 14, 30}, 'power', 1, ...
           'beam', {3, 4, 5}, 'rIdx', {100, 200, 150}, 'dIdx', {10, 20, 30}, ...
           'pulse', {'short', 'long', 'long'}, 'hits', {80, 60, 70});
ppi_update(ppi, sys, 100, 12.0, 2, p);
pos = getpixelposition(ppi.ax);
xl0 = get(ppi.ax, 'XLim');
tol = 8 * max(diff(xl0)/pos(3), diff(get(ppi.ax, 'YLim'))/pos(4));
txtNear = ppi_callbacks('hover', fig, plot_xy(p(3)) + [0.75*tol 0]);
xy = @(k) plot_xy(p(k));
S = getappdata(fig, 'ppi');
ok = true;
check = @(cond, msg) report(cond, msg);

%% Data tips
txt = ppi_callbacks('hover', fig, xy(1) + [0.05 0.05]);
ok = check(~isempty(strfind(txt, 'short pulse')) && ~isempty(strfind(txt, '10.000 km')) ...
           && ~isempty(strfind(txt, '12.50 m/s')) && ~isempty(strfind(txt, 'no DBF')), ...
           'hover over the short-pulse plot shows its values') && ok;
ok = check(strcmp(get(S.h.tip, 'Visible'), 'on'), 'data tip is visible') && ok;
txt = ppi_callbacks('hover', fig, [-60 -60]);
ok = check(isempty(txt) && strcmp(get(S.h.tip, 'Visible'), 'off'), 'no data tip far from all plots') && ok;
tgtXY = (cfg.targets.pos0_m(1:2) + cfg.targets.vel_mps(1:2)*12) / 1e3;
txt = ppi_callbacks('hover', fig, tgtXY);
ok = check(~isempty(strfind(txt, 'True target T1')), 'hover over a true target') && ok;

ok = check(~isempty(strfind(txtNear, '30.000 km')), 'data tip within 8 px at the full view') && ok;

%% Zoom (wheel and keys) and pan
xl0 = get(ppi.ax, 'XLim');
ppi_callbacks('zoom', fig, 2, xy(3));
xl = get(ppi.ax, 'XLim');
c = xy(3);
ok = check(abs(diff(xl) - diff(xl0)/2) < 1e-9 && abs((c(1) - xl(1))/diff(xl) - (c(1) - xl0(1))/diff(xl0)) < 1e-9, ...
           'zoom x2 around a point keeps the point in place') && ok;
txt = ppi_callbacks('hover', fig, xy(3) + [0.75*tol 0]);   % 6 px at full view, 12 px zoomed
ok = check(isempty(txt), 'data tip tolerance (8 px) follows the zoom') && ok;
ppi_callbacks('pan', fig, [5 -3]);
ok = check(all(abs(get(ppi.ax, 'XLim') - (xl + 5)) < 1e-9), 'pan moves the view') && ok;
ppi_callbacks('key', fig, struct('Key', 'add', 'Character', '+'));
ok = check(abs(diff(get(ppi.ax, 'XLim')) - diff(xl)/1.5) < 1e-9, 'key + zooms in') && ok;
ppi_callbacks('key', fig, struct('Key', 'r', 'Character', 'r'));
ok = check(all(abs(get(ppi.ax, 'XLim') - xl0) < 1e-9), 'key r resets the view') && ok;
ppi_callbacks('zoom', fig, 1e6, [0 0]);
ok = check(abs(diff(xl0)/diff(get(ppi.ax, 'XLim')) - 500) < 1e-6, 'zoom is limited to x500') && ok;
ppi_callbacks('zoom', fig, 1e-6, [0 0]);
ok = check(all(abs(get(ppi.ax, 'XLim') - xl0) < 1e-6), 'zoom out stops at the full view') && ok;

%% Filters
ppi_callbacks('option', fig, 'showShort', false);
ok = check(isempty(get(S.h.short, 'XData')) && isempty(ppi_callbacks('hover', fig, xy(1))), ...
           'hiding the short-pulse plots removes them and their data tips') && ok;
ppi_callbacks('option', fig, 'showShort', true);
ppi_callbacks('option', fig, 'minSnr_dB', 20);
ok = check(numel(get(S.h.long, 'XData')) == 1, 'SNR filter hides the 14 dB plot') && ok;
ppi_callbacks('option', fig, 'minSnr_dB', 0);
ppi_callbacks('option', fig, 'history', 1);                % last scan: 10 s at 6 rpm
ok = check(numel(get(S.h.long, 'XData')) == 1, 'age filter hides the plot of the previous scan') && ok;
ppi_callbacks('option', fig, 'history', 0);
ok = check(numel(get(S.h.long, 'XData')) == 2, 'all plots shown again') && ok;
ppi_callbacks('option', fig, 'showLabels', true);
S = getappdata(fig, 'ppi');
ok = check(numel(S.h.labels) == 3, 'one label per plot') && ok;
ppi_callbacks('option', fig, 'showLabels', false);
ppi_callbacks('option', fig, 'showRings', false);
ok = check(all(strcmp(get(S.h.rings, 'Visible'), 'off')), 'range rings hidden') && ok;
ppi_callbacks('option', fig, 'showRings', true);

%% Pinned data tips
n = ppi_callbacks('pin', fig, xy(2));
ok = check(n == 1, 'click on a plot pins its data tip') && ok;
n = ppi_callbacks('pin', fig, xy(3));
ok = check(n == 2, 'second pinned data tip') && ok;
n = ppi_callbacks('pin', fig, xy(2));
ok = check(n == 1, 'click on a pinned plot removes its data tip') && ok;
ppi_callbacks('cleartips', fig);
S = getappdata(fig, 'ppi');
ok = check(isempty(S.h.pins), 'remove all data tips') && ok;

%% Pause, clear, save
ppi_callbacks('pause', fig);
ppi_update(ppi, sys, 120, 13.0, 2, p(1));
ok = check(numel(get(S.h.long, 'XData')) == 2, 'paused display is not redrawn') && ok;
ppi_callbacks('pause', fig);
ok = check(isempty(get(S.h.long, 'XData')), 'resumed display shows the new data') && ok;
ppi_update(ppi, sys, 120, 13.0, 2, p);
ppi_callbacks('clear', fig);
ok = check(isempty(get(S.h.short, 'XData')) && isempty(get(S.h.long, 'XData')), 'clear removes the plots') && ok;
file = [tempname() '.png'];
ppi_callbacks('save', fig, file);
ok = check(exist(file, 'file') == 2, 'PPI saved as PNG') && ok;
if exist(file, 'file'), delete(file); end
close(fig);
fprintf('Interactive PPI: %s\n', ternary(ok, 'PASS', 'FAIL'));
end

function v = plot_xy(p)
% Ground position of a plot on the PPI [km].
el = p.el_deg;
if isnan(el), el = 0; end
v = p.range_m * cosd(el) * [sind(p.az_deg) cosd(p.az_deg)] / 1e3;
end

function ok = report(cond, msg)
ok = logical(cond);
fprintf('  %-62s %s\n', msg, ternary(ok, 'ok', 'FAIL'));
end

function s = ternary(c, a, b)
if c, s = a; else, s = b; end
end
