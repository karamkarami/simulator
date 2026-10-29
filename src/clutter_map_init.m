function cm = clutter_map_init(sys, wv)
%CLUTTER_MAP_INIT  Empty clutter map of sub-pulse wv (cfg.clutterRejection.method = 'map').
%   cm = CLUTTER_MAP_INIT(sys, wv) returns the map of the zero-Doppler band
%   (sys.cmap.bins) of every azimuth cell, with nothing learned yet:
%     M     [nBeams x nRange x nBins x nAz] mean power of every cell
%     count [1 x nAz] scans averaged into every azimuth cell (0 = not observed)
%     cur, sum, n: the azimuth cell of the current pass of the antenna and the
%     sum of its n band powers (they enter the map when the pass is over, see
%     clutter_map_fold.m)
nBins = numel(sys.cmap.bins);
cm.M     = zeros(sys.nBeams, wv.nRange, nBins, sys.cmap.nAz, 'single');
cm.count = zeros(1, sys.cmap.nAz);
cm.cur   = 0;
cm.sum   = zeros(sys.nBeams, wv.nRange, nBins);
cm.n     = 0;
end
