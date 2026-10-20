function adcData = adc_sample(x, sys)
%ADC_SAMPLE  ADC model: ADC noise (sets the effective number of bits), rounding and clipping.
%   adcData = ADC_SAMPLE(x, sys) converts the analog signal x (in LSB) to int16.
if sys.cfg.adc.addAdcNoise
    x = x + sys.adc.sigmaAdc * randn(size(x), 'single');
end
adcData = int16(x);                                   % rounds and clips to 16 bits
if sys.cfg.adc.nBits < 16
    fs = sys.adc.fullScale;
    adcData = min(max(adcData, -fs), fs - 1);
end
end
