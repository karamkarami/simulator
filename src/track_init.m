function tr = track_init(params)
%TRACK_INIT  Empty EKF tracker (see track_update.m).
tr.p = params;
tr.tracks = struct('id', {}, 'confirmed', {}, 'x', {}, 'P', {}, 't', {}, 'hits', {}, 'misses', {});
tr.nextId = 1;
tr.log = zeros(0, 11);   % one row per update: [time id confirmed x y z vx vy vz hits plotIndex]
tr.nPlots = 0;
end
