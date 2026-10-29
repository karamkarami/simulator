function cm = clutter_map_fold(cm, sys)
%CLUTTER_MAP_FOLD  The pass of the antenna over an azimuth cell enters the clutter map.
%   cm = CLUTTER_MAP_FOLD(cm, sys) adds the mean band power of the CPI
%   updates of the current pass (cm.sum / cm.n) to the map of azimuth cell
%   cm.cur by exponential averaging, map = map + w*(power - map), with
%   w = max(sys.cmap.weight, 1/n) for the n-th scan of the cell (the plain
%   mean of the first scans). Called when the antenna enters another azimuth
%   cell, at a skipped sector and at the end of the run.
if cm.cur > 0 && cm.n > 0
    a = cm.cur;
    n = cm.count(a) + 1;
    w = max(sys.cmap.weight, 1/n);
    cm.M(:, :, :, a) = cm.M(:, :, :, a) + w * (single(cm.sum / cm.n) - cm.M(:, :, :, a));
    cm.count(a) = min(n, numel(sys.cmap.alpha));          % (selects the threshold factor)
end
cm.cur = 0;
cm.sum(:) = 0;
cm.n = 0;
end
