function [isDet, noise] = cfar_2d(P, sys)
%CFAR_2D  Two-dimensional cell-averaging CFAR on every range-Doppler map.
%   [isDet, noise] = CFAR_2D(P, sys) tests every cell of P [nBeams x nRange x nDoppler].
%   The noise of a cell is the mean of the training cells of a rectangular
%   window (range x Doppler) around it, excluding the guard cells. The Doppler
%   axis is circular. Cells closer than the window size to the first or last
%   range cell, and the Doppler bins within +-notchBins of zero (MTI notch,
%   where the clutter residue is), are not tested.

c  = sys.cfar;
oR = c.guardRange + c.trainRange;      % half size of the outer window in range
oD = c.guardDoppler + c.trainDoppler;  % half size of the outer window in Doppler
sumOuter = box_sum(P, oR, oD);
sumGuard = box_sum(P, c.guardRange, c.guardDoppler);
noise = (sumOuter - sumGuard) / c.nTrain;
isDet = P > c.alpha * noise;
[~, nR, nD] = size(P);
isDet(:, [1:oR, nR-oR+1:nR], :) = false;
zeroBin = floor(nD/2) + 1;                                     % zero Doppler after fftshift
isDet(:, :, zeroBin-c.notchBins:zeroBin+c.notchBins) = false;
end

function S = box_sum(P, hr, hd)
% Sum of P over a (2*hr+1) x (2*hd+1) window centred on every cell (range x Doppler).
S = convn(P, ones(1, 2*hr+1), 'same');                         % range (zero padded)
S = cat(3, S(:, :, end-hd+1:end), S, S(:, :, 1:hd));           % circular padding in Doppler
S = convn(S, ones(1, 1, 2*hd+1), 'valid');
end
