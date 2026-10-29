function [isDet, cm] = clutter_map_detect(cm, isDet, P, az_deg, sys)
%CLUTTER_MAP_DETECT  Clutter-map test of the zero-Doppler band of one CPI update.
%   [isDet, cm] = CLUTTER_MAP_DETECT(cm, isDet, P, az_deg, sys) takes the CFAR
%   decisions isDet and the power map P [nBeams x nRange x nDoppler] of a
%   CPI update whose newest pulse has the antenna azimuth az_deg. A cell of
%   the band sys.cmap.bins stays a detection only if its power also exceeds
%   alpha(n) times the map of its azimuth cell, n = the scans averaged into
%   it (sys.cmap.alpha); cells of an azimuth cell that no earlier scan has
%   observed are not detections. The band powers are added to the current
%   pass over the azimuth cell; when the antenna enters another cell, the
%   pass enters the map (clutter_map_fold.m). The map is thus built from the
%   earlier scans only, so a target does not raise its own threshold.
%   The azimuth is rounded to single precision as in the UDP header, so the
%   C receiver puts every update into the same azimuth cell.
x = mod(double(single(az_deg)), 360);
a = min(floor(x / sys.cmap.azCell_deg), sys.cmap.nAz - 1) + 1;    % azimuth cell of this update
if a ~= cm.cur
    cm = clutter_map_fold(cm, sys);                       % the antenna has left the previous cell
    cm.cur = a;
end
band = sys.cmap.bins;
Pb = P(:, :, band);
n = cm.count(a);
if n == 0
    isDet(:, :, band) = false;                            % this cell has no map yet
else
    isDet(:, :, band) = isDet(:, :, band) & Pb > sys.cmap.alpha(n) * cm.M(:, :, :, a);
end
cm.sum = cm.sum + double(Pb);
cm.n = cm.n + 1;
end
