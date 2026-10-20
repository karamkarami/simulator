function y = fir_decimate_valid(x, h, D)
%FIR_DECIMATE_VALID  FIR filtering with a valid convolution followed by decimation.
%   y = FIR_DECIMATE_VALID(x, h, D) filters every column of x with h and keeps
%   every D-th sample of the valid part of the convolution. It gives the same
%   result as
%       y = conv(x, h, 'valid');   y = y(1:D:end, :);
%   but only computes the samples that are kept (about D times faster).

[N, nCol] = size(x);
L = numel(h);
g = flipud(h(:));                                   % y(k) = sum_i g(i) * x(start_k + i)
K = floor((N - L)/D) + 1;                           % number of output samples
P = ceil(L/D);                                      % taps per polyphase branch
g(P*D) = 0;                                         % zero-pad the taps to P*D
G = reshape(g, D, P).';                             % row p: taps (p-1)*D+1 ... p*D
M = K + P - 1;                                      % blocks of D input samples that are used
if M*D > N
    x(M*D, :) = 0;                                  % only multiplied by zero-padded taps
end
X = reshape(x(1:M*D, :), D, M*nCol);                % column m: samples (m-1)*D+1 ... m*D
Z = reshape(G * X, P, M, nCol);                     % every branch applied to every block
y = reshape(Z(1, 1:K, :), K, nCol);
for p = 2:P
    y = y + reshape(Z(p, p:p+K-1, :), K, nCol);     % output k uses blocks k ... k+P-1
end
end
