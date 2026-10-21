function out = ppi_callbacks(action, fig, varargin)
%PPI_CALLBACKS  Interaction with the PPI created by ppi_init.m.
%   out = PPI_CALLBACKS(action, fig, ...) is called by the figure and its
%   controls; it can also be called directly (tests, scripts):
%     'redraw'                    draw the stored data with the current options
%     'scroll', evt               mouse wheel: zoom around the cursor
%     'zoom', factor, [x y]       zoom by factor (> 1 = in) around [x y] km
%     'pan', [dx dy]              move the view by [dx dy] km
%     'down' / 'move' / 'up'      mouse buttons and motion (drag = pan, click = pin)
%     'hover', [x y]              data tip of the object nearest to [x y] km
%                                 (out = its text, '' if nothing is close)
%     'pin', [x y]                pin (or remove) the data tip of the nearest object
%     'key', evt                  + - arrows r c space
%     'option', name, value       showShort, showLong, showTruth, showRings,
%                                 showLabels, showAntenna, minSnr_dB, history
%     'reset' 'clear' 'cleartips' 'pause'   buttons of the panel
%     'save' [, fileName]         save the PPI as a PNG file
out = [];
if ~ishandle(fig), return; end
S = getappdata(fig, 'ppi');
switch action
    case 'redraw'
        S = redraw(S);
    case 'scroll'
        evt = varargin{1};
        p = cursor_km(S);
        if ~isempty(p)
            S = zoom_at(S, 1.25^(-evt.VerticalScrollCount), p);
        end
    case 'zoom'
        S = zoom_at(S, varargin{1}, varargin{2});
    case 'pan'
        S = set_limits(S, get(S.ax, 'XLim') + varargin{1}(1), get(S.ax, 'YLim') + varargin{1}(2));
    case 'down'
        p = cursor_km(S);
        if isempty(p), return; end                          % outside the axes
        if strcmp(get(fig, 'SelectionType'), 'open')        % double click
            S = set_limits(S, S.home(1:2), S.home(3:4));
            S.drag = [];
        else
            S.drag = struct('px', figure_px(fig), 'xl', get(S.ax, 'XLim'), 'yl', get(S.ax, 'YLim'), ...
                            'moved', false, 'p', p);
        end
    case 'move'
        if ~isempty(S.drag)                                 % pan with the mouse
            d = figure_px(fig) - S.drag.px;
            if any(abs(d) > 3), S.drag.moved = true; end
            if S.drag.moved
                k = km_per_pixel(S, S.drag.xl, S.drag.yl);
                S = set_limits(S, S.drag.xl - d(1)*k, S.drag.yl - d(2)*k);
            end
        else
            p = cursor_km(S);
            S = show_cursor(S, p);
            [S, out] = hover(S, p);
        end
    case 'up'
        if ~isempty(S.drag) && ~S.drag.moved                % a click: pin the data tip
            S = pin(S, S.drag.p);
        end
        S.drag = [];
    case 'hover'
        [S, out] = hover(S, varargin{1});
    case 'pin'
        S = pin(S, varargin{1});
        out = numel(S.h.pins);
    case 'key'
        S = key(S, varargin{1});
    case 'option'
        S.opt.(varargin{1}) = varargin{2};
        S = redraw(S);
    case 'reset'
        S = set_limits(S, S.home(1:2), S.home(3:4));
    case 'clear'
        S.opt.clearTime_s = S.data.t_s;
        S = redraw(S);
    case 'cleartips'
        S = clear_pins(S);
    case 'pause'
        S.opt.paused = ~S.opt.paused;
        labels = {'Pause display', 'Resume display'};
        set(S.h.ui.pause, 'String', labels{1 + S.opt.paused});
        if ~S.opt.paused, S = redraw(S); end
    case 'save'
        if ~isempty(varargin)
            file = varargin{1};
        else
            [name, folder] = uiputfile('*.png', 'Save the PPI', 'ppi.png');
            if isequal(name, 0), return; end
            file = fullfile(folder, name);
        end
        setappdata(fig, 'ppi', S);
        print(fig, file, '-dpng');
        out = file;
        return
    otherwise
        error('ppi_callbacks: unknown action ''%s''', action);
end
setappdata(fig, 'ppi', S);
end

%% ------------------------------------------------------------------------
function S = redraw(S)
% Antenna, true targets and the plots that pass the filters.
h = S.h;  o = S.opt;  d = S.data;  sys = S.sys;
onOff = {'off', 'on'};
set(h.ant, 'XData', [0 S.rMax*sind(d.az_deg)], 'YData', [0 S.rMax*cosd(d.az_deg)], ...
           'Visible', onOff{1 + o.showAntenna});
set(h.rings, 'Visible', onOff{1 + o.showRings});

targets = sys.cfg.targets;
if o.showTruth && ~isempty(targets)
    n = numel(targets);
    xy = zeros(n, 2);
    for i = 1:n
        p = targets(i).pos0_m + targets(i).vel_mps*d.t_s;
        xy(i, :) = p(1:2) / 1e3;
        p0 = targets(i).pos0_m / 1e3;
        set(h.names(i), 'Position', [xy(i, :) 0], 'Visible', 'on');
        set(h.tracks(i), 'XData', [p0(1) xy(i, 1)], 'YData', [p0(2) xy(i, 2)], 'Visible', 'on');
    end
    set(h.truth, 'XData', xy(:, 1), 'YData', xy(:, 2), 'Visible', 'on');
else
    set(h.truth, 'Visible', 'off');
    set(h.names, 'Visible', 'off');
    set(h.tracks, 'Visible', 'off');
end

% Plots that pass the filters
plots = d.plots;
nAll = numel(plots);
x = zeros(1, nAll);  y = x;  show = false(1, nAll);  isShort = show;
if nAll > 0
    isShort = strcmp({plots.pulse}, 'short');
    el = [plots.el_deg];
    el(isnan(el)) = 0;                                     % no elevation without beamforming
    groundRange = [plots.range_m] .* cosd(el) / 1e3;
    x = groundRange .* sind([plots.az_deg]);
    y = groundRange .* cosd([plots.az_deg]);
    t = [plots.time_s];
    show = ((isShort & o.showShort) | (~isShort & o.showLong)) ...
         & [plots.snr_dB] >= o.minSnr_dB & t > o.clearTime_s;
    if o.history > 0
        show = show & (d.t_s - t) < o.history * 360 / sys.ant.rate_dps;
    end
end
s = show & isShort;  l = show & ~isShort;
set(h.short, 'XData', x(s), 'YData', y(s));
set(h.long,  'XData', x(l), 'YData', y(l));
set_datatip_rows(h.short, plots(s));
set_datatip_rows(h.long,  plots(l));
S.vis = struct('x', x(show), 'y', y(show), 'idx', find(show));

% Labels next to the plots
delete(h.labels(ishandle(h.labels)));
S.h.labels = [];
if o.showLabels
    idx = find(show);
    for k = idx(1:min(end, 300))
        S.h.labels(end+1) = text(x(k), y(k), sprintf('  %.1f km, %+.1f m/s', plots(k).range_m/1e3, ...
            plots(k).vr_mps), 'Parent', S.ax, 'Color', [0.8 0.8 0.8], 'FontSize', 8, 'HitTest', 'off', 'Clipping', 'on');
    end
end
set(h.title, 'String', sprintf('Scan %d   t = %.2f s   az = %5.1f deg   plots %d / %d   zoom x%.1f', ...
    d.scan, d.t_s, d.az_deg, sum(show), nAll, zoom_factor(S)));
end

function set_datatip_rows(hLine, plots)
% MATLAB data tips (axes toolbar) with the plot values (not in Octave).
if isempty(plots) || ~isprop(hLine, 'DataTipTemplate'), return; end
try
    tpl = hLine.DataTipTemplate;
    el = [plots.el_deg];
    tpl.DataTipRows = [dataTipTextRow('Range [km]', [plots.range_m]/1e3, '%.3f'), ...
                       dataTipTextRow('Azimuth [deg]', [plots.az_deg], '%.2f'), ...
                       dataTipTextRow('Elevation [deg]', el, '%.2f'), ...
                       dataTipTextRow('Radial vel. [m/s]', [plots.vr_mps], '%.2f'), ...
                       dataTipTextRow('SNR [dB]', [plots.snr_dB], '%.1f'), ...
                       dataTipTextRow('Time [s]', [plots.time_s], '%.3f'), ...
                       dataTipTextRow('CPI hits', [plots.hits], '%d')];
catch
end
end

%% ------------------------------------------------------------------------
function [S, txt] = hover(S, p)
% Data tip of the plot or true target nearest to p (within 8 pixels).
txt = '';
[txt, pos] = nearest_text(S, p);
if isempty(txt)
    set(S.h.tip, 'Visible', 'off');
    set(S.h.sel, 'XData', NaN, 'YData', NaN);
    return
end
k = km_per_pixel(S, get(S.ax, 'XLim'), get(S.ax, 'YLim'));
xl = get(S.ax, 'XLim');
if pos(1) > mean(xl)                                       % keep the tip inside the axes
    set(S.h.tip, 'Position', [pos(1) - 10*k, pos(2) + 6*k, 0], 'HorizontalAlignment', 'right');
else
    set(S.h.tip, 'Position', [pos(1) + 10*k, pos(2) + 6*k, 0], 'HorizontalAlignment', 'left');
end
set(S.h.tip, 'String', txt, 'Visible', 'on');
set(S.h.sel, 'XData', pos(1), 'YData', pos(2));
end

function [txt, pos] = nearest_text(S, p)
% Text and position of the nearest plot or true target, '' if none is close.
txt = '';  pos = [];
if isempty(p), return; end
tol = 8 * km_per_pixel(S, get(S.ax, 'XLim'), get(S.ax, 'YLim'));
best = tol;
for j = 1:numel(S.vis.x)                                  % plots
    dist = hypot(S.vis.x(j) - p(1), S.vis.y(j) - p(2));
    if dist <= best
        best = dist;
        pos = [S.vis.x(j) S.vis.y(j)];
        txt = plot_text(S.data.plots(S.vis.idx(j)));
    end
end
if S.opt.showTruth                                         % true targets
    xt = get(S.h.truth, 'XData');  yt = get(S.h.truth, 'YData');
    for i = 1:numel(xt)
        dist = hypot(xt(i) - p(1), yt(i) - p(2));
        if dist < best
            best = dist;
            pos = [xt(i) yt(i)];
            txt = truth_text(S, i);
        end
    end
end
end

function txt = plot_text(p)
if isnan(p.el_deg), el = '- (no DBF)'; else, el = sprintf('%.2f deg', p.el_deg); end
txt = sprintf(['Plot (%s pulse)\nRange       %.3f km\nAzimuth     %.2f deg\nElevation   %s\n' ...
               'Radial vel. %.2f m/s\nSNR         %.1f dB\nTime        %.3f s\nCPI hits    %d\n' ...
               'Beam %d, range cell %d, Doppler bin %d'], p.pulse, p.range_m/1e3, p.az_deg, el, ...
               p.vr_mps, p.snr_dB, p.time_s, p.hits, p.beam, p.rIdx, p.dIdx);
end

function txt = truth_text(S, i)
tgt = S.sys.cfg.targets(i);
[R, az, el, vr] = target_geometry(tgt, S.data.t_s, S.sys.ant.height_m);
pos = tgt.pos0_m + tgt.vel_mps * S.data.t_s;
txt = sprintf(['True target %s\nRange       %.3f km\nAzimuth     %.2f deg\nElevation   %.2f deg\n' ...
               'Height      %.0f m\nRadial vel. %.1f m/s\nSpeed       %.1f m/s\nMean RCS    %.2f m^2'], ...
               tgt.name, R/1e3, az, el, pos(3), vr, norm(tgt.vel_mps), tgt.rcs_m2);
end

function S = pin(S, p)
% Pin the data tip of the nearest object, or remove a pinned tip at p.
k = km_per_pixel(S, get(S.ax, 'XLim'), get(S.ax, 'YLim'));
for i = 1:numel(S.h.pins)                                  % click on a pinned tip: remove it
    if hypot(S.h.pins{i}.x - p(1), S.h.pins{i}.y - p(2)) <= 8*k
        delete([S.h.pins{i}.txt S.h.pins{i}.mk]);
        S.h.pins(i) = [];
        return
    end
end
[txt, pos] = nearest_text(S, p);
if isempty(txt), return; end
pinTxt = text(pos(1) + 10*k, pos(2) - 6*k, txt, 'Parent', S.ax, 'Color', 'w', ...
              'BackgroundColor', [0.05 0.2 0.05], 'EdgeColor', [0.3 0.9 0.3], 'FontSize', 9, ...
              'VerticalAlignment', 'top', 'Margin', 4, 'HitTest', 'off');
mk = plot(S.ax, pos(1), pos(2), 's', 'Color', [0.3 0.9 0.3], 'MarkerSize', 13, 'LineWidth', 1.5, 'HitTest', 'off', 'Clipping', 'on');
S.h.pins{end+1} = struct('txt', pinTxt, 'mk', mk, 'x', pos(1), 'y', pos(2));
end

function S = clear_pins(S)
for i = 1:numel(S.h.pins)
    delete([S.h.pins{i}.txt S.h.pins{i}.mk]);
end
S.h.pins = {};
end

function S = key(S, evt)
xl = get(S.ax, 'XLim');  yl = get(S.ax, 'YLim');
step = 0.1 * diff(xl);
c = [mean(xl) mean(yl)];
switch evt.Key
    case {'add', 'equal', 'plus'},       S = zoom_at(S, 1.5, c);
    case {'subtract', 'hyphen', 'minus'}, S = zoom_at(S, 1/1.5, c);
    case 'leftarrow',  S = set_limits(S, xl - step, yl);
    case 'rightarrow', S = set_limits(S, xl + step, yl);
    case 'uparrow',    S = set_limits(S, xl, yl + step);
    case 'downarrow',  S = set_limits(S, xl, yl - step);
    case {'r', 'home'}, S = set_limits(S, S.home(1:2), S.home(3:4));
    case 'c',           S = clear_pins(S);
    case 'space'
        S.opt.paused = ~S.opt.paused;
        labels = {'Pause display', 'Resume display'};
        set(S.h.ui.pause, 'String', labels{1 + S.opt.paused});
        if ~S.opt.paused, S = redraw(S); end
end
end

%% ------------------------------------------------------------------------
function S = zoom_at(S, factor, c)
% Zoom by factor around the point c [km] (the point stays under the cursor).
xl = get(S.ax, 'XLim');  yl = get(S.ax, 'YLim');
span = diff(S.home(1:2));
f = max(min(factor, diff(xl)/(span/500)), diff(xl)/span); % between full view and x500
S = set_limits(S, c(1) + (xl - c(1))/f, c(2) + (yl - c(2))/f);
end

function S = set_limits(S, xl, yl)
set(S.ax, 'XLim', xl, 'YLim', yl);
set(S.h.title, 'String', regexprep(get(S.h.title, 'String'), 'zoom x[\d.]+', ...
    sprintf('zoom x%.1f', zoom_factor(S))));
end

function z = zoom_factor(S)
z = diff(S.home(1:2)) / diff(get(S.ax, 'XLim'));
end

function k = km_per_pixel(S, xl, yl)
% Data units per screen pixel (axis equal: the larger of the two scales).
pos = getpixelposition(S.ax);
k = max(diff(xl)/pos(3), diff(yl)/pos(4));
end

function p = cursor_km(S)
% Cursor position in the axes [km], [] when it is outside the axes.
cp = get(S.ax, 'CurrentPoint');
p = cp(1, 1:2);
xl = get(S.ax, 'XLim');  yl = get(S.ax, 'YLim');
if p(1) < xl(1) || p(1) > xl(2) || p(2) < yl(1) || p(2) > yl(2), p = []; end
end

function S = show_cursor(S, p)
if isempty(p)
    set(S.h.cursor, 'String', '');
    return
end
az = mod(atan2(p(1), p(2)) * 180/pi, 360);
set(S.h.cursor, 'String', {sprintf('Cursor: %.2f km, %.1f deg', hypot(p(1), p(2)), az), ...
                           sprintf('E %.2f km, N %.2f km', p(1), p(2))});
end

function px = figure_px(fig)
% Mouse position in the figure in pixels.
units = get(fig, 'Units');
set(fig, 'Units', 'pixels');
px = get(fig, 'CurrentPoint');
set(fig, 'Units', units);
px = px(1:2);
end
