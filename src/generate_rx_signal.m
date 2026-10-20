function x = generate_rx_signal(sys, targets, t_s, azAnt_deg, clutterIf)
%GENERATE_RX_SIGNAL  Received IF signal of all elements during the receive record of one PRI.
%   x = GENERATE_RX_SIGNAL(sys, targets, t_s, azAnt_deg, clutterIf) returns
%   the real signal [nAdc x nElements] at the ADC input, in ADC units (LSB),
%   before quantization: target echoes of both sub-pulses on the 30 MHz IF +
%   clutter + thermal noise, set to zero while the long pulse is transmitted
%   (receiver blanked). t_s is the start time of the PRI and azAnt_deg the
%   antenna azimuth. clutterIf is the output of clutter_signal.m ([] = no clutter).

c   = sys.cfg.const.c;
fs  = sys.fs;
% Single precision is plenty for a 16-bit ADC and about twice as fast
x   = zeros(sys.nAdc, sys.nCh, 'single');
names = {'short', 'long'};

%% Target echoes
for i = 1:numel(targets)
    [R, az, el, vr] = target_geometry(targets(i), t_s, sys.ant.height_m);
    azRel = wrap180(az - azAnt_deg);
    if abs(azRel) > sys.ant.azSpan_deg
        continue                                   % far outside the azimuth beam
    end
    [v, u] = steering_vector(sys, azRel, el);
    elRel  = asind(u);
    % Two-way amplitude factor: product of the one-way power patterns
    gain = pattern_gain(sys.ant.az_deg, sys.ant.azGain_dB, azRel)^2 ...
         * pattern_gain(sys.ant.elRel_deg, sys.ant.txElGain_dB, elRel) ...
         * pattern_gain(sys.ant.elRel_deg, sys.ant.elElemGain_dB, elRel);
    power = sys.radarConst * targets(i).rcs_m2 / R^4 * gain;   % per element [W]
    amp   = sqrt(2*power) * sys.adc.scale;                    % IF amplitude [LSB]
    for p = 1:2
        wv   = sys.wave.(names{p});
        tau0 = 2*R/c;
        % ADC samples that contain the echo
        first = max(1, ceil((wv.tTx_s + tau0 - sys.rxStart_s)*fs) + 1);
        last  = min(sys.nAdc, floor((wv.tTx_s + tau0 + wv.width_s - sys.rxStart_s)*fs) + 1);
        if first > last, continue; end
        n  = (first:last)';
        t  = sys.tAdc(n);                           % time inside the PRI
        tau = 2*(R + vr*t)/c;                       % delay of every sample (intra-pulse Doppler)
        % RF echo mixed down by the LO: IF carrier plus the carrier phase -2*pi*fLO*tau
        echo = amp * waveform_envelope(wv, t - wv.tTx_s - tau0) ...
             .* exp(1j*2*pi*(sys.fIF*(t - wv.tTx_s - tau) - sys.fLO*tau));
        x(n, :) = x(n, :) + real(single(echo) * single(v).');
    end
end

%% Clutter (same echo on all elements, different phase)
if ~isempty(clutterIf)
    x = x + real(single(clutterIf) * single(sys.clutterSteer).');
end

%% Thermal noise
if sys.cfg.rx.addThermalNoise
    x = x + sys.adc.sigmaTh * randn(sys.nAdc, sys.nCh, 'single');
end

%% Receiver blanked while the long pulse is transmitted
x = x .* sys.rxGate;
end
