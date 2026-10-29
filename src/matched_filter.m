function rp = matched_filter(beams, wv)
%MATCHED_FILTER  Pulse compression of every beam (valid convolution).
%   rp = MATCHED_FILTER(beams, wv) filters every row of beams [nBeams x nSamples]
%   (the DDC samples wv.ddcIdx) with the pulse compression filter of sub-pulse
%   wv (matched or mismatch) and returns the range profiles rp [nBeams x nRange]
%   on the range axis wv.range_m (delay axis wv.delay_s, index axis
%   1:wv.nRange). The samples outside the receive interval count as zeros
%   (wv.pad, see setup_radar.m).
nB = size(beams, 1);
rp = conv2([zeros(nB, wv.pad(1)), beams, zeros(nB, wv.pad(2))], wv.mfTaps, 'valid');
end
