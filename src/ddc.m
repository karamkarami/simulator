function bb = ddc(adcData, sys, wv)
%DDC  Two-stage digital down-converter.
%   bb = DDC(adcData, sys, wv) converts the real ADC samples of all channels
%   ([nAdc x nElements]) to complex baseband for sub-pulse wv:
%     stage 1: mix the IF wv.fIF_Hz to 0 Hz, low-pass filter h1, decimate by 10
%     stage 2: low-pass filter h2, decimate by 6
%   Both filters use valid convolutions (no zero padding at the edges). Only
%   the first wv.nAdcUsed samples, needed for the processed range interval of
%   this sub-pulse, are converted. Output: [nSamples x nElements].
x  = single(adcData(1:wv.nAdcUsed, :));               % single precision: enough for 16-bit data
bb = fir_decimate_valid(x,  sys.ddc.h1, sys.ddc.decim(1), wv.fIF_Hz / sys.fs);
bb = fir_decimate_valid(bb, sys.ddc.h2, sys.ddc.decim(2));
end
