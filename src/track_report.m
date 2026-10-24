function rep = track_report(tracker)
%TRACK_REPORT  The current tracks for the display (one struct per track).
%   Fields: id, confirmed, t_s (last update), pos_m [x y z] and vel_mps
%   [vx vy vz] relative to the antenna (x East, y North, z up; z = 0 in 2-D),
%   alt_m (height above ground, NaN in 2-D), hits, misses, dim.
tk = tracker.tracks;
d = tracker.p.dim;
rep = struct('id', {}, 'confirmed', {}, 't_s', {}, 'pos_m', {}, 'vel_mps', {}, 'alt_m', {}, ...
             'hits', {}, 'misses', {}, 'dim', {});
for i = 1:numel(tk)
    pos = zeros(1, 3);  vel = zeros(1, 3);
    pos(1:d) = tk(i).x(1:d);
    vel(1:d) = tk(i).x(d+1:end);
    alt = NaN;
    if d == 3, alt = pos(3) + tracker.p.height_m; end
    rep(i) = struct('id', tk(i).id, 'confirmed', tk(i).confirmed, 't_s', tk(i).t, 'pos_m', pos, ...
                    'vel_mps', vel, 'alt_m', alt, 'hits', tk(i).hits, 'misses', tk(i).misses, 'dim', d);
end
end
