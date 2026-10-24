function tr = track_update(tr, plots, t_s)
%TRACK_UPDATE  Extended Kalman filter tracking of the plots.
%   tr = TRACK_UPDATE(tr, plots, t_s) processes the new plots (in their
%   order) and then deletes the tracks that were not updated for too long at
%   time t_s (call it every PRI, also without plots). tr comes from
%   track_init(sys.track).
%
%   State relative to the antenna (x East, y North, z up): 2-D [x y vx vy],
%   3-D [x y z vx vy vz]; nearly constant velocity: between two updates the
%   state is predicted with F = [I dt*I; 0 I] and the white-acceleration
%   process noise q*[dt^3/3 dt^2/2; dt^2/2 dt] per axis.
%   Measurement of a plot: 2-D z = [ground range; azimuth] (ground range =
%   range*cos(elevation) when the elevation is known), 3-D z = [range;
%   azimuth; elevation], with the nonlinear h(x) linearised at the predicted
%   state (EKF; the azimuth innovation is wrapped to +-pi).
%   Association: the plot goes to the track with the smallest d^2 + ln|S|
%   among the tracks whose gate (d^2 = nu'*inv(S)*nu <= p.gate) holds it; a
%   plot without a track starts a tentative track at the measured position
%   with zero velocity and velocity std p.sigmaV0. A tentative track with
%   p.confirmHits updates is confirmed. Missed scans: the antenna passed the
%   predicted position half a scan ago without an update; a track with more
%   than p.maxMiss (tentative: p.maxMissTentative) missed scans is deleted.
%   The updates are appended to tr.log ([time id confirmed x y z vx vy vz hits
%   plotIndex], z = 0 in 2-D).

p = tr.p;
n = 2 * p.dim;
for k = 1:numel(plots)
    tr.nPlots = tr.nPlots + 1;
    [z, ok] = measurement(plots(k), p.dim);
    if ~ok, continue; end
    tp = plots(k).time_s;
    best = 0;  bestCost = inf;
    for j = 1:numel(tr.tracks)
        [x, P] = predict(tr.tracks(j).x, tr.tracks(j).P, tp - tr.tracks(j).t, p);
        [nu, H] = innovation(z, x, p.dim);
        S = H*P*H' + p.R;
        d2 = nu' * (S \ nu);
        if d2 <= p.gate
            cost = d2 + log(det(S));
            if cost < bestCost
                best = j;  bestCost = cost;
                xb = x;  Pb = P;  nub = nu;  Hb = H;  Sb = S;
            end
        end
    end
    if best > 0                                         % EKF update (Joseph form)
        K = Pb * Hb' / Sb;
        I_KH = eye(n) - K*Hb;
        tk = tr.tracks(best);
        tk.x = xb + K*nub;
        tk.P = I_KH*Pb*I_KH' + K*p.R*K';
        tk.P = (tk.P + tk.P') / 2;
        tk.t = tp;
        tk.hits = tk.hits + 1;
        tk.misses = 0;
        tk.confirmed = tk.confirmed || tk.hits >= p.confirmHits;
        tr.tracks(best) = tk;
    else                                                % new tentative track
        [pos, J] = position(z, p.dim);
        P = zeros(n);
        P(1:p.dim, 1:p.dim) = J*p.R*J';
        P(p.dim+1:end, p.dim+1:end) = p.sigmaV0^2 * eye(p.dim);
        tk = struct('id', tr.nextId, 'confirmed', p.confirmHits <= 1, 'x', [pos; zeros(p.dim, 1)], ...
                    'P', P, 't', tp, 'hits', 1, 'misses', 0);
        tr.nextId = tr.nextId + 1;
        tr.tracks(end+1) = tk;
        best = numel(tr.tracks);
    end
    tk = tr.tracks(best);
    xyz = zeros(1, 6);
    xyz([1:p.dim, 4:3+p.dim]) = tk.x.';
    tr.log(end+1, :) = [tk.t, tk.id, tk.confirmed, xyz, tk.hits, tr.nPlots];
end

% Missed scans and deletion
if isempty(tr.tracks), return; end
miss = max(0, floor((t_s - [tr.tracks.t]) / p.scanTime_s - 0.5));
for j = 1:numel(tr.tracks)
    tr.tracks(j).misses = miss(j);
end
limit = p.maxMissTentative * ones(size(miss));
limit([tr.tracks.confirmed]) = p.maxMiss;
tr.tracks(miss > limit) = [];
end

function [z, ok] = measurement(pl, dim)
% Measurement vector of a plot.
el = deg2rad(pl.el_deg);
az = deg2rad(pl.az_deg);
if dim == 3
    z = [pl.range_m; az; el];
    ok = isfinite(el);
else
    rho = pl.range_m;
    if isfinite(el), rho = rho * cos(el); end
    z = [rho; az];
    ok = true;
end
end

function [x, P] = predict(x, P, dt, p)
% Nearly-constant-velocity prediction over dt (may be slightly negative).
d = p.dim;
F = [eye(d), dt*eye(d); zeros(d), eye(d)];
a = abs(dt);
Q = [diag(p.q)*a^3/3, diag(p.q)*a^2/2; diag(p.q)*a^2/2, diag(p.q)*a];
x = F*x;
P = F*P*F' + Q;
end

function [nu, H] = innovation(z, x, dim)
% Innovation z - h(x) and Jacobian H of h at x.
px = x(1);  py = x(2);
rh2 = px^2 + py^2;  rh = sqrt(rh2);
if dim == 3
    pz = x(3);
    R2 = rh2 + pz^2;  R = sqrt(R2);
    h = [R; atan2(px, py); atan2(pz, rh)];
    H = [px/R,          py/R,          pz/R,    0 0 0;
         py/rh2,        -px/rh2,       0,       0 0 0;
         -px*pz/(R2*rh), -py*pz/(R2*rh), rh/R2, 0 0 0];
else
    h = [rh; atan2(px, py)];
    H = [px/rh,  py/rh,  0 0;
         py/rh2, -px/rh2, 0 0];
end
nu = z - h;
nu(2) = mod(nu(2) + pi, 2*pi) - pi;                 % azimuth wrap
end

function [pos, J] = position(z, dim)
% Position of a measurement and the Jacobian d(position)/d(z).
az = z(2);
if dim == 3
    R = z(1);  el = z(3);
    pos = R * [cos(el)*sin(az); cos(el)*cos(az); sin(el)];
    J = [cos(el)*sin(az),  R*cos(el)*cos(az), -R*sin(el)*sin(az);
         cos(el)*cos(az), -R*cos(el)*sin(az), -R*sin(el)*cos(az);
         sin(el),          0,                  R*cos(el)];
else
    rho = z(1);
    pos = rho * [sin(az); cos(az)];
    J = [sin(az),  rho*cos(az);
         cos(az), -rho*sin(az)];
end
end
