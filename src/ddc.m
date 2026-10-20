function bb = ddc(adcData, sys)
%DDC  Two-stage digital down-converter of all channels.
%   bb = DDC(adcData, sys) converts the real ADC samples [nAdc x nElements]
%   on the 30 MHz IF to complex baseband [nDdc x nElements]:
%     mixer   : real line x .* cos(2*pi*fIF*t), imaginary line x .* (-sin(2*pi*fIF*t))
%     stage 1 : low-pass filter h1 and decimation by 10 (120 -> 12 MHz)
%     stage 2 : low-pass filter h2 and decimation by 6  (12 -> 2 MHz)
%   Both filters use valid convolutions (no zero padding at the edges). The
%   time of every output sample is sys.tDdc; the samples of the short and of
%   the long sub-pulse are sys.wave.short.ddcIdx and sys.wave.long.ddcIdx.
x  = single(adcData);                              % single precision: enough for 16-bit data
xI = x .* sys.ddc.cos;                             % real (I) line
xQ = x .* sys.ddc.minusSin;                        % imaginary (Q) line
h1 = sys.ddc.h1;  D1 = sys.ddc.decim(1);
h2 = sys.ddc.h2;  D2 = sys.ddc.decim(2);
I  = fir_decimate_valid(fir_decimate_valid(xI, h1, D1), h2, D2);
Q  = fir_decimate_valid(fir_decimate_valid(xQ, h1, D1), h2, D2);
bb = complex(I, Q);
end
