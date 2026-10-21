function ppi = ppi_init(sys)
%PPI_INIT  Create the interactive PPI (plan position indicator) display.
%   ppi = PPI_INIT(sys) creates the figure. ppi_update.m draws the antenna,
%   the true targets and the plots; ppi_callbacks.m handles the interaction:
%     mouse wheel          zoom in/out around the cursor
%     left drag            pan
%     double click         reset the view (whole instrumented range)
%     hover over a plot    data tip with range, azimuth, elevation, radial
%                          velocity, SNR, time, sub-pulse, hits, beam
%     click on a plot      pin its data tip (click again to remove it)
%     keys                 + / - zoom, arrows pan, r reset, c clear the
%                          pinned data tips, space pause/resume the display
%   The panel on the right shows or hides the short/long-pulse plots, the
%   true targets, range rings, plot labels and the antenna; filters the plots
%   by SNR and age; clears the plots; and saves the PPI as a PNG file.
%   The MATLAB zoom/pan/data tip tools of the axes toolbar also work.
cfg  = sys.cfg;
rMax = cfg.display.maxRange_m;
if isempty(rMax), rMax = cfg.rx.maxRange_m; end
rMax = rMax / 1e3;                                         % km
if cfg.display.visible, vis = 'on'; else, vis = 'off'; end

fig = figure('Name', 'PPI', 'Color', 'k', 'Visible', vis, 'InvertHardcopy', 'off', ...
             'NumberTitle', 'off', 'Position', [100 100 1000 760]);
ax  = axes('Parent', fig, 'Color', 'k', 'XColor', [0.6 0.6 0.6], 'YColor', [0.6 0.6 0.6], ...
           'Units', 'normalized', 'Position', [0.07 0.12 0.66 0.8]);
hold(ax, 'on');
axis(ax, 'equal');
axis(ax, rMax*[-1.05 1.05 -1.05 1.05]);
xlabel(ax, 'East [km]');
ylabel(ax, 'North [km]');
if exist('disableDefaultInteractivity', 'file')
    disableDefaultInteractivity(ax);                       % our own wheel zoom / drag pan
end

%% Range rings and azimuth spokes
ringStep = 10^floor(log10(rMax/2));
if rMax/ringStep > 6, ringStep = 2*ringStep; end
phi  = linspace(0, 2*pi, 361);
gray = [0.25 0.25 0.25];
h.rings = [];
for r = ringStep:ringStep:rMax
    h.rings(end+1) = plot(ax, r*sin(phi), r*cos(phi), '-', 'Color', gray, 'HitTest', 'off', 'Clipping', 'on');
    h.rings(end+1) = text(r*sind(45), r*cosd(45), sprintf('%g km', r), 'Color', [0.5 0.5 0.5], ...
                          'Parent', ax, 'HitTest', 'off', 'Clipping', 'on');
end
for a = 0:30:330
    h.rings(end+1) = plot(ax, [0 rMax*sind(a)], [0 rMax*cosd(a)], '-', 'Color', gray, 'HitTest', 'off', 'Clipping', 'on');
    h.rings(end+1) = text(1.03*rMax*sind(a), 1.03*rMax*cosd(a), sprintf('%d', a), 'Color', [0.5 0.5 0.5], ...
                          'Parent', ax, 'HorizontalAlignment', 'center', 'HitTest', 'off', 'Clipping', 'on');
end

%% Antenna, true targets, plots
h.ant    = plot(ax, [0 0], [0 rMax], '-', 'Color', [0 0.8 0], 'LineWidth', 1.5, 'HitTest', 'off', 'Clipping', 'on');
h.truth  = plot(ax, NaN, NaN, 's', 'Color', [1 0.3 0.3], 'MarkerSize', 10);
h.names  = gobjects_compat(numel(cfg.targets));
h.tracks = gobjects_compat(numel(cfg.targets));
for i = 1:numel(cfg.targets)
    h.names(i)  = text(NaN, NaN, ['  ' cfg.targets(i).name], 'Color', [1 0.3 0.3], 'Parent', ax, 'HitTest', 'off', 'Clipping', 'on');
    h.tracks(i) = plot(ax, NaN, NaN, ':', 'Color', [1 0.3 0.3], 'HitTest', 'off', 'Clipping', 'on');   % true track since t = 0
end
h.short = plot(ax, NaN, NaN, 'o', 'Color', [0 0.9 1], 'MarkerFaceColor', [0 0.9 1], 'MarkerSize', 5);
h.long  = plot(ax, NaN, NaN, 'o', 'Color', [1 0.9 0], 'MarkerFaceColor', [1 0.9 0], 'MarkerSize', 5);
h.labels = [];                                             % plot labels (option)
hLeg = legend([h.truth h.short h.long], {'true target', 'plot (short pulse)', 'plot (long pulse)'}, ...
       'TextColor', 'w', 'Color', 'k', 'Location', 'southoutside', 'Orientation', 'horizontal');
try
    set(hLeg, 'AutoUpdate', 'off');                        % data tips and labels stay out of the legend
catch
end
h.title = title(ax, '', 'Color', 'w');

%% Data tip (hover) and selection marker
h.tip = text(0, 0, '', 'Parent', ax, 'Color', 'w', 'BackgroundColor', [0.12 0.12 0.12], ...
             'EdgeColor', [0.8 0.8 0.8], 'FontSize', 9, 'VerticalAlignment', 'bottom', ...
             'Visible', 'off', 'HitTest', 'off', 'Margin', 4);
h.sel = plot(ax, NaN, NaN, 'o', 'Color', 'w', 'MarkerSize', 13, 'LineWidth', 1.5, 'HitTest', 'off', 'Clipping', 'on');
h.pins = {};                                               % pinned data tips: {text, marker}

%% Control panel (right)
bg = [0.15 0.15 0.15];  fg = [0.9 0.9 0.9];
x0 = 0.76;  w = 0.22;  rowH = 0.04;
y = 0.92;
uicontrol(fig, 'Style', 'text', 'String', 'Display', 'Units', 'normalized', 'Position', [x0 y w rowH], ...
          'BackgroundColor', 'k', 'ForegroundColor', 'w', 'FontWeight', 'bold', 'HorizontalAlignment', 'left');
checks = {'showShort', 'Short-pulse plots', true;  'showLong', 'Long-pulse plots', true;
          'showTruth', 'True targets', logical(cfg.display.showTruth);  'showRings', 'Range rings', true;
          'showLabels', 'Plot labels', false;  'showAntenna', 'Antenna', true};
for i = 1:size(checks, 1)
    y = y - rowH;
    h.ui.(checks{i, 1}) = uicontrol(fig, 'Style', 'checkbox', 'String', checks{i, 2}, 'Value', checks{i, 3}, ...
        'Units', 'normalized', 'Position', [x0 y w rowH], 'BackgroundColor', bg, 'ForegroundColor', fg, ...
        'Callback', @(src, evt) ppi_callbacks('option', fig, checks{i, 1}, logical(get(src, 'Value'))));
    opt.(checks{i, 1}) = checks{i, 3};
end
y = y - 1.5*rowH;
uicontrol(fig, 'Style', 'text', 'String', 'Min. SNR [dB]', 'Units', 'normalized', 'Position', [x0 y 0.12 rowH], ...
          'BackgroundColor', 'k', 'ForegroundColor', fg, 'HorizontalAlignment', 'left');
h.ui.minSnr = uicontrol(fig, 'Style', 'edit', 'String', '0', 'Units', 'normalized', ...
    'Position', [x0+0.13 y 0.09 rowH], 'BackgroundColor', bg, 'ForegroundColor', fg, ...
    'Callback', @(src, evt) ppi_callbacks('option', fig, 'minSnr_dB', str2double(get(src, 'String'))));
opt.minSnr_dB = 0;
y = y - rowH;
uicontrol(fig, 'Style', 'text', 'String', 'Plots shown', 'Units', 'normalized', 'Position', [x0 y 0.12 rowH], ...
          'BackgroundColor', 'k', 'ForegroundColor', fg, 'HorizontalAlignment', 'left');
h.ui.history = uicontrol(fig, 'Style', 'popupmenu', 'String', {'all', 'last scan', 'last 2 scans'}, ...
    'Value', 1, 'Units', 'normalized', 'Position', [x0+0.13 y 0.09 rowH], 'BackgroundColor', bg, ...
    'ForegroundColor', fg, 'Callback', @(src, evt) ppi_callbacks('option', fig, 'history', get(src, 'Value') - 1));
opt.history = 0;                                           % 0 all, 1 last scan, 2 last two scans
buttons = {'Reset view', 'reset';  'Clear plots', 'clear';  'Remove data tips', 'cleartips';
           'Pause display', 'pause';  'Save PNG ...', 'save'};
y = y - 0.5*rowH;
for i = 1:size(buttons, 1)
    y = y - 1.2*rowH;
    h.ui.(buttons{i, 2}) = uicontrol(fig, 'Style', 'pushbutton', 'String', buttons{i, 1}, 'Units', 'normalized', ...
        'Position', [x0 y w rowH], 'BackgroundColor', bg, 'ForegroundColor', fg, ...
        'Callback', @(src, evt) ppi_callbacks(buttons{i, 2}, fig));
end
y = y - 2.2*rowH;
h.cursor = uicontrol(fig, 'Style', 'text', 'String', '', 'Units', 'normalized', 'Position', [x0 y w 2*rowH], ...
                     'BackgroundColor', 'k', 'ForegroundColor', [0.7 0.9 0.7], 'HorizontalAlignment', 'left');
uicontrol(fig, 'Style', 'text', 'Units', 'normalized', 'Position', [x0 0.01 w 0.2], ...
    'BackgroundColor', 'k', 'ForegroundColor', [0.6 0.6 0.6], 'HorizontalAlignment', 'left', ...
    'String', {'Wheel: zoom', 'Drag: pan', 'Double click: reset', 'Hover a plot: data tip', ...
               'Click a plot: pin the tip', 'Keys: + - arrows r c space'});

%% State (used by ppi_update and ppi_callbacks)
opt.clearTime_s = -Inf;                                    % plots up to this time are cleared
opt.paused = false;
S.h = h;
S.ax = ax;
S.rMax = rMax;
S.home = rMax*[-1.05 1.05 -1.05 1.05];
S.opt = opt;
S.sys = sys;
S.data = struct('az_deg', 0, 't_s', 0, 'scan', 1, 'plots', []);
S.vis = struct('x', [], 'y', [], 'idx', []);               % visible plots (for picking)
S.drag = [];
setappdata(fig, 'ppi', S);

set(fig, 'WindowScrollWheelFcn',  @(src, evt) ppi_callbacks('scroll', fig, evt), ...
         'WindowButtonDownFcn',   @(src, evt) ppi_callbacks('down', fig), ...
         'WindowButtonMotionFcn', @(src, evt) ppi_callbacks('move', fig), ...
         'WindowButtonUpFcn',     @(src, evt) ppi_callbacks('up', fig), ...
         'KeyPressFcn',           @(src, evt) ppi_callbacks('key', fig, evt));

ppi.fig = fig;
ppi.ax = ax;
ppi.rMax = rMax;
end

function g = gobjects_compat(n)
% Empty graphics handle array (gobjects in MATLAB, doubles in Octave).
if exist('gobjects', 'builtin') || exist('gobjects', 'file')
    g = gobjects(1, n);
else
    g = zeros(1, n);
end
end
