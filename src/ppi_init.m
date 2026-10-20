function ppi = ppi_init(sys)
%PPI_INIT  Create the PPI (plan position indicator) display.
cfg  = sys.cfg;
rMax = cfg.display.maxRange_m;
if isempty(rMax), rMax = cfg.rx.maxRange_m; end
rMax = rMax / 1e3;                                         % km
if cfg.display.visible, vis = 'on'; else, vis = 'off'; end

ppi.fig = figure('Name', 'PPI', 'Color', 'k', 'Visible', vis, 'InvertHardcopy', 'off');
ppi.ax  = axes('Parent', ppi.fig, 'Color', 'k', 'XColor', [0.6 0.6 0.6], 'YColor', [0.6 0.6 0.6]);
hold(ppi.ax, 'on');
axis(ppi.ax, 'equal');
axis(ppi.ax, rMax*[-1.05 1.05 -1.05 1.05]);
xlabel(ppi.ax, 'East [km]');
ylabel(ppi.ax, 'North [km]');

% Range rings and azimuth spokes
ringStep = 10^floor(log10(rMax/2));
if rMax/ringStep > 6, ringStep = 2*ringStep; end
phi = linspace(0, 2*pi, 361);
grid = [0.25 0.25 0.25];
for r = ringStep:ringStep:rMax
    plot(ppi.ax, r*sin(phi), r*cos(phi), '-', 'Color', grid);
    text(r*sind(45), r*cosd(45), sprintf('%g km', r), 'Color', [0.5 0.5 0.5], 'Parent', ppi.ax);
end
for a = 0:30:330
    plot(ppi.ax, [0 rMax*sind(a)], [0 rMax*cosd(a)], '-', 'Color', grid);
end

ppi.rMax   = rMax;
ppi.hAnt   = plot(ppi.ax, [0 0], [0 rMax], '-', 'Color', [0 0.8 0], 'LineWidth', 1.5);
ppi.hTruth = plot(ppi.ax, NaN, NaN, 's', 'Color', [1 0.3 0.3], 'MarkerSize', 10);
ppi.hNames  = cell(1, numel(cfg.targets));
ppi.hTracks = cell(1, numel(cfg.targets));
for i = 1:numel(cfg.targets)
    ppi.hNames{i}  = text(NaN, NaN, ['  ' cfg.targets(i).name], 'Color', [1 0.3 0.3], 'Parent', ppi.ax);
    ppi.hTracks{i} = plot(ppi.ax, NaN, NaN, ':', 'Color', [1 0.3 0.3]);     % true track since t = 0
end
ppi.hShort = plot(ppi.ax, NaN, NaN, 'o', 'Color', [0 0.9 1], 'MarkerFaceColor', [0 0.9 1], 'MarkerSize', 5);
ppi.hLong  = plot(ppi.ax, NaN, NaN, 'o', 'Color', [1 0.9 0], 'MarkerFaceColor', [1 0.9 0], 'MarkerSize', 5);
legend([ppi.hTruth ppi.hShort ppi.hLong], {'true target', 'plot (short pulse)', 'plot (long pulse)'}, ...
       'TextColor', 'w', 'Color', 'k', 'Location', 'southoutside', 'Orientation', 'horizontal');
ppi.title = title(ppi.ax, '', 'Color', 'w');
end
