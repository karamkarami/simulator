function rcs = target_rcs(tgt, scan)
%TARGET_RCS  RCS of a target during a scan.
%   With Swerling 1 (see draw_target_rcs.m) every scan has its own RCS,
%   constant during the scan; otherwise the mean RCS tgt.rcs_m2 is used.
if isfield(tgt, 'rcsScan') && ~isempty(tgt.rcsScan)
    rcs = tgt.rcsScan(min(scan, numel(tgt.rcsScan)));
else
    rcs = tgt.rcs_m2;
end
end
