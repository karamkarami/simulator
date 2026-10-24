function p = track_setup(cfg, sys)
%TRACK_SETUP  Parameters of the EKF tracker (track_update.m) from cfg.track.
%   p.dim is 3 (range, azimuth, elevation measured; needs the elevation of
%   the DBF) or 2 (range and azimuth in the horizontal plane). With
%   cfg.track.dim = 'auto' it is 3 with beamforming and 2 without.
t = cfg.track;
p.enable = logical(t.enable);
switch lower(t.dim)
    case 'auto', p.dim = 2 + sys.dbf.enabled;
    case '2d',   p.dim = 2;
    case '3d'
        p.dim = 3;
        if ~sys.dbf.enabled
            warning('cfg.track.dim = ''3d'' needs the elevation of the DBF (1 element: no DBF): tracking in 2-D.');
            p.dim = 2;
        end
    otherwise
        error('cfg.track.dim must be ''auto'', ''2d'' or ''3d''.');
end
sig = [t.sigmaRange_m, deg2rad(t.sigmaAz_deg), deg2rad(t.sigmaEl_deg)];
p.R = diag(sig(1:p.dim).^2);                       % measurement noise [range az (el)]
q = [t.sigmaAccel_mps2, t.sigmaAccel_mps2, t.sigmaAccelZ_mps2].^2;
p.q = q(1:p.dim);                                  % white-acceleration spectral densities per axis
p.sigmaV0 = t.maxSpeed_mps / sqrt(3);              % velocity std of a new track (uniform in +-maxSpeed)
p.gate = 2 * gammaincinv(t.gateProb, p.dim/2);     % chi-square gate with dim degrees of freedom
p.confirmHits = t.confirmHits;
p.maxMissTentative = t.maxMissTentative;
p.maxMiss = t.maxMiss;
p.scanTime_s = 360 / sys.ant.rate_dps;
p.height_m = sys.ant.height_m;                     % the state is relative to the antenna
end
