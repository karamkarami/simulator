function ppi_update(ppi, sys, az_deg, t_s, scan, plots)
%PPI_UPDATE  Refresh the PPI: antenna position, true targets and extracted plots.
%   The data is stored in the figure, so the interactive options of
%   ppi_callbacks.m (filters, zoom, data tips) can redraw it at any time.
%   While the display is paused (space key or the Pause button) the data is
%   still stored but not drawn.
if ~ishandle(ppi.fig), return; end                       % the window was closed
S = getappdata(ppi.fig, 'ppi');
S.sys = sys;
S.data = struct('az_deg', az_deg, 't_s', t_s, 'scan', scan, 'plots', plots);
setappdata(ppi.fig, 'ppi', S);
if ~S.opt.paused
    ppi_callbacks('redraw', ppi.fig);
end
drawnow;
end
