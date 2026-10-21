function adcData = adc_sample(x, sys)
%ADC_SAMPLE  16-bit ADC with 12.7 effective bits.
%   adcData = ADC_SAMPLE(x, sys) converts the ADC input x (sqrt(W), see
%   generate_rx_signal.m) to LSB, adds the ADC noise that sets the effective
%   number of bits, rounds and clips to int16. With cfg.adc.quantize = false
%   (SNR tests only) the analog values in LSB are returned instead.
x = x * sys.adc.lsbPerUnit;                           % [sqrt(W)] -> LSB
if sys.cfg.adc.addAdcNoise
    x = x + sys.adc.sigmaAdc * randn(size(x), 'single');
end
if ~sys.cfg.adc.quantize
    adcData = x;
    return
end
adcData = int16(x);                                   % rounds and clips to 16 bits
if sys.cfg.adc.nBits < 16
    fs = sys.adc.fullScale;
    adcData = min(max(adcData, -fs), fs - 1);
end
end
