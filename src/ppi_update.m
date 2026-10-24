function sys = ppi_update(ppi, sys, az_deg, t_s, scan, plots, dets)
%PPI_UPDATE  Refresh the PPI: antenna position, true targets, plots and detections.
%   sys = PPI_UPDATE(ppi, sys, az_deg, t_s, scan, plots, dets) stores the
%   data in the figure, so the interactive options of ppi_callbacks.m
%   (filters, zoom, data tips) can redraw it at any time. dets are the CFAR
%   detections before the plot extractor (rows of detection_log_init.m;
%   optional). While the display is paused (space key or the Pause button)
%   the data is still stored but not drawn.
%   When the user applied a new CFAR threshold in the panel, the returned sys
%   holds it (sys.cfar.alpha), so the processing uses it from now on.
if nargin < 7, dets = zeros(0, 10); end
if ~ishandle(ppi.fig), return; end                       % the window was closed
dB = ppi_callbacks('takeCfar', ppi.fig);
if ~isempty(dB)
    sys.cfar.alpha = 10^(dB/10);
    if sys.cfg.sim.verbose
        fprintf('  CFAR threshold set to %.2f dB (Pfa %.1e) at t = %.3f s\n', dB, ...
                cfar_pfa(sys.cfar.alpha, sys.cfar.nTrain, sys.dop.nLooks), t_s);
    end
end
S = getappdata(ppi.fig, 'ppi');
S.sys = sys;
S.data = struct('az_deg', az_deg, 't_s', t_s, 'scan', scan, 'plots', plots, 'dets', dets);
setappdata(ppi.fig, 'ppi', S);
if ~isempty(dB)
    ppi_callbacks('cfarApplied', ppi.fig, dB);
end
if ~S.opt.paused
    ppi_callbacks('redraw', ppi.fig);
end
drawnow;
end
