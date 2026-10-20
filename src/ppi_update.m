function ppi_update(ppi, sys, az_deg, t_s, scan, plots)
%PPI_UPDATE  Refresh the PPI: antenna position, true targets and extracted plots.
set(ppi.hAnt, 'XData', [0 ppi.rMax*sind(az_deg)], 'YData', [0 ppi.rMax*cosd(az_deg)]);

if sys.cfg.display.showTruth && ~isempty(sys.cfg.targets)
    n = numel(sys.cfg.targets);
    xy = zeros(n, 2);
    for i = 1:n
        p = sys.cfg.targets(i).pos0_m + sys.cfg.targets(i).vel_mps*t_s;
        xy(i, :) = p(1:2) / 1e3;
    end
    set(ppi.hTruth, 'XData', xy(:, 1), 'YData', xy(:, 2));
    for i = 1:n
        p0 = sys.cfg.targets(i).pos0_m / 1e3;
        set(ppi.hNames{i}, 'Position', [xy(i, :) 0]);
        set(ppi.hTracks{i}, 'XData', [p0(1) xy(i, 1)], 'YData', [p0(2) xy(i, 2)]);
    end
end

if ~isempty(plots)
    isShort = strcmp({plots.pulse}, 'short');
    groundRange = [plots.range_m] .* cosd([plots.el_deg]) / 1e3;
    x = groundRange .* sind([plots.az_deg]);
    y = groundRange .* cosd([plots.az_deg]);
    set(ppi.hShort, 'XData', x(isShort),  'YData', y(isShort));
    set(ppi.hLong,  'XData', x(~isShort), 'YData', y(~isShort));
end
set(ppi.title, 'String', sprintf('Scan %d   t = %.2f s   az = %5.1f deg   plots = %d', ...
    scan, t_s, az_deg, numel(plots)));
drawnow;
end
