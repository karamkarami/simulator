function y = fir_decimate_valid(x, h, D)
%FIR_DECIMATE_VALID  FIR filter (valid time-domain convolution) followed by decimation.
%   y = FIR_DECIMATE_VALID(x, h, D) slides the filter h over every column of
%   x and keeps every D-th output of the valid part of the convolution:
%       y(k) = sum_j h(j) * x((k-1)*D + L + 1 - j),   j = 1 ... L = numel(h)
%   which is the same as conv(x, h, 'valid') followed by y(1:D:end), but only
%   the kept outputs are computed. The sum is done tap by tap: for tap j all
%   output positions are updated at once (no FFT).

[N, nCol] = size(x);
L = numel(h);
K = floor((N - L)/D) + 1;                           % number of output samples
y = zeros(K, nCol);
for j = 1:L
    y = y + h(L + 1 - j) * x(j : D : j + (K-1)*D, :);
end
end
