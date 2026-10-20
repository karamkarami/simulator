function u = waveform_envelope(wv, t)
%WAVEFORM_ENVELOPE  Complex envelope of a sub-pulse.
%   u = WAVEFORM_ENVELOPE(wv, t) evaluates the envelope at the times t [s]
%   (t = 0 is the start of the pulse). The envelope is zero outside the pulse.
%   wv.type = 'lfm'       : linear FM chirp from -B/2 to +B/2
%   wv.type = 'phasecode' : P4 polyphase code with chips of wv.chip_s seconds

switch lower(wv.type)
    case 'lfm'
        k = wv.bandwidth_Hz / wv.width_s;               % chirp rate [Hz/s]
        u = exp(1j*pi*k*(t - wv.width_s/2).^2);
    case 'phasecode'
        nChips = round(wv.width_s / wv.chip_s);
        m = min(max(floor(t / wv.chip_s), 0), nChips-1); % chip index (0-based)
        u = exp(1j*(pi*m.^2/nChips - pi*m));             % P4 code phases
    otherwise
        error('Unknown waveform type "%s"', wv.type);
end
u(t < 0 | t >= wv.width_s) = 0;
end
