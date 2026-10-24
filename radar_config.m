function cfg = radar_config()
%RADAR_CONFIG  All user-editable parameters of the rotating search radar simulator.
%
%   cfg = RADAR_CONFIG() returns a struct holding every parameter of the
%   simulation. Edit the values below, or change fields of cfg after calling
%   this function (this is what the tests do), e.g.
%
%       cfg = radar_config();
%       cfg.prf.mode = 'staggered';
%       cfg.clutter.enable = false;
%       results = run_radar_sim(cfg);
%
%   Quantities that come from input .mat files (antenna patterns, DDC
%   filters, DBF matrix, staggered PRIs and MTI weights) are listed under
%   "Input files". Default files are created by make_input_files.m when they
%   are missing.

root = fileparts(mfilename('fullpath'));
addpath(fullfile(root, 'src'));

%% Physical constants
cfg.const.c  = 299792458;       % speed of light [m/s]
cfg.const.kB = 1.380649e-23;    % Boltzmann constant [J/K]
cfg.const.T0 = 290;             % reference noise temperature [K]

%% Input files
cfg.files.antenna = fullfile(root, 'inputs', 'antenna_pattern.mat'); % patterns, element spacing, gains
cfg.files.ddc     = fullfile(root, 'inputs', 'ddc_filters.mat');     % low-pass filters of the 2 DDC stages
cfg.files.dbf     = fullfile(root, 'inputs', 'dbf_matrix.mat');      % elevation beamforming matrix
cfg.files.stagger = fullfile(root, 'inputs', 'stagger_prf.mat');     % staggered PRIs + MTI weights

%% Radar site and antenna
cfg.fc            = 3e9;  % center frequency [Hz]
cfg.ant.nElements = 48;   % receive elements in elevation (1 = single antenna, beamforming disabled)
cfg.ant.elementSpacing_lambda = 0.66;  % element spacing in elevation [wavelengths]
% (azimuth: fixed 4.1 deg beam of every element; patterns come from the antenna file)
cfg.ant.height_m  = 6;    % antenna height above the ground [m]
cfg.ant.tilt_deg  = 11;   % array tilt above the horizon [deg]
cfg.ant.rpm       = 6;    % rotation speed [rev/min]
cfg.ant.az0_deg   = 0;    % antenna azimuth at t = 0, clockwise from North [deg]

%% Transmitter and receiver
cfg.tx.peakPower_W     = 20e3;   % peak transmit power [W]
cfg.rx.noiseFigure_dB  = 4;      % receiver noise figure [dB]
cfg.rx.losses_dB       = 3;      % system losses not modelled elsewhere (RF, radome, ...) [dB]
cfg.rx.gain_dB         = 60;     % gain of the receiver path from the antenna to the ADC [dB]
cfg.rx.txIsolation_dB  = 40;     % transmitter-to-receiver isolation (leakage while transmitting) [dB]
cfg.rx.maxRange_m      = 80e3;   % instrumented range of the long pulse: sets the receive window [m]
cfg.rx.addThermalNoise = true;   % false only for debugging / SNR measurements
cfg.rx.addLeakage      = true;   % false only for debugging / SNR measurements

%% ADC
cfg.adc.fs_Hz            = 120e6;  % sampling frequency [Hz]
cfg.adc.nBits            = 16;     % ADC word length [bit]
cfg.adc.enob             = 12.7;   % effective number of bits
cfg.adc.fullScale_dBm    = 10;     % power of a full-scale sine at the ADC input (2 Vpp into 50 ohm) [dBm]
cfg.adc.addAdcNoise      = true;   % false only for debugging / SNR measurements
cfg.adc.quantize         = true;   % false only for SNR measurements (keeps the analog values)
cfg.adc.fIF_Hz           = 30e6;   % IF carrier of the received signal (receiver LO = fc - 30 MHz)

%% Waveforms: every PRI holds a short and a long sub-pulse (short first)
% PRI timeline: short pulse at its tTx -> receive (near range) -> long pulse
% at its tTx -> receive (far range) until the end of the record. While a
% pulse is transmitted the receiver sees the transmitter leakage. Both
% sub-pulses are on the same 30 MHz IF and are separated in time. The long
% pulse is blind up to c*150us/2 = 22.5 km; with tTx = 200 us the short
% pulse covers 4.2 ... 23.7 km and hands over to the long pulse.
cfg.wave.short.type         = 'lfm';    % 'lfm' or 'phasecode' (P4 polyphase code)
cfg.wave.short.width_s      = 25e-6;    % pulse width [s]
cfg.wave.short.bandwidth_Hz = 1e6;      % LFM sweep bandwidth [Hz]
cfg.wave.short.chip_s       = 1e-6;     % phase-code chip length [s]
cfg.wave.short.tTx_s        = 0;        % transmit time inside the PRI [s]

cfg.wave.long.type          = 'lfm';
cfg.wave.long.width_s       = 150e-6;
cfg.wave.long.bandwidth_Hz  = 1e6;
cfg.wave.long.chip_s        = 1e-6;
cfg.wave.long.tTx_s         = 200e-6;

%% PRF
cfg.prf.mode  = 'fixed';   % 'fixed' or 'staggered' (staggered PRIs come from cfg.files.stagger)
cfg.prf.pri_s = 1e-3;      % PRI in fixed mode [s]

%% DDC (filters come from cfg.files.ddc)
cfg.ddc.decimation = [10 6];   % decimation factor of stage 1 and stage 2 (120 -> 12 -> 2 MHz)
cfg.ddc.coefFracBits = 13;     % filter coefficients in the file: signed 16-bit, 13 fractional bits
cfg.ddc.shift      = [13 13];  % right shift of the accumulator of each stage (13 = unit gain);
                               % every stage output is rounded and saturated to signed 16 bit
cfg.ddc.fixedPoint = true;     % false only for SNR measurements (no rounding / saturation)

%% Matched filter and Doppler processing
cfg.proc.taylorNbar   = 4;     % Taylor window (LFM only): number of nearly constant sidelobes
cfg.proc.taylorSll_dB = -35;   % Taylor window peak sidelobe level [dB]
cfg.proc.cpiLength    = [];    % pulses per CPI buffer ([] = pulses arriving while a target is in the
                               % 3 dB azimuth beam: round(beamwidth / (rotation rate * PRI)) = 114)
cfg.proc.nfft         = [];    % Doppler FFT length ([] = next power of 2: of the pulses in the buffer
                               % with a fixed PRF (114 -> 128), of the samples per row when staggered (28 -> 32))
cfg.proc.updateEvery  = 1;     % run MTI/Doppler/CFAR every N new pulses (1 = every pulse)

%% 2-D CA-CFAR on every range-Doppler map
cfg.cfar.pfa          = 1e-8;  % design probability of false alarm (per cell and CPI)
cfg.cfar.guardRange   = 3;     % guard cells on each side in range
cfg.cfar.trainRange   = 8;     % training cells on each side in range
cfg.cfar.guardDoppler = 3;     % guard cells on each side in Doppler
cfg.cfar.trainDoppler = 4;     % training cells on each side in Doppler
cfg.cfar.notchWidth   = 1;     % half-width of the zero-Doppler band that is not tested (MTI notch,
                               % clutter residue), in Doppler window mainlobes (1 -> +-2 bins)

%% Pruning of the CFAR detections
cfg.prune.sidelobeMargin_dB = 20;  % a peak this much weaker than a peak on the same range, beam
                                   % or Doppler line is treated as its sidelobe and removed

%% Plot extraction (azimuth peak picking over successive CPIs)
cfg.extract.minExtent_deg = 2; % minimum azimuth extent of the detections of a plot [deg]
                               % (noise false alarms last only a few CPIs, targets about one beamwidth)
cfg.extract.maxMiss = 2;       % CPIs without a detection before a plot is closed

%% Ground clutter
cfg.clutter.enable           = true;
cfg.clutter.sigma0_dB        = -30;   % land clutter reflectivity [dB m^2/m^2]
cfg.clutter.minRange_m       = 1000;  % clutter starts here [m]
cfg.clutter.maxRange_m       = [];    % [] = radar horizon (4/3 earth model)
cfg.clutter.azStep_deg       = 0.25;  % azimuth size of a clutter cell [deg]
cfg.clutter.rangeStepSamples = 12;    % range size of a clutter cell in ADC samples (12 -> 15 m)
cfg.clutter.icmStd_mps       = 0.1;   % std of internal clutter motion (wind) [m/s]
cfg.clutter.elevation_deg    = 0;     % elevation of the clutter returns [deg]

%% Simulation length
cfg.sim.nScans      = 1;      % number of antenna revolutions
cfg.sim.sectors_deg = [];     % [] = full 360 deg; else rows [azStart azEnd]: only these sectors are simulated
cfg.sim.seed        = 1;      % random seed
cfg.sim.swerling    = 1;      % target RCS model: 1 = Swerling 1 (new RCS every scan), 0 = constant
cfg.sim.verbose     = true;

%% Targets (see make_target.m)
cfg.targets = six_targets_scenario();

%% Recording of the raw ADC data and of the DDC output, streaming of the DDC output
cfg.io.outFolder = fullfile(root, 'output');
cfg.io.saveRaw   = true;                                    % raw ADC data: one .mat file per PRI
cfg.io.rawFolder = fullfile(root, 'output', 'raw_pulses');  % -> pulse_<n>.mat (adcData, hdr)
cfg.io.saveDdc   = false;                                   % DDC output: one .mat file per PRI
cfg.io.ddcFolder = fullfile(root, 'output', 'ddc_pulses');  % -> ddc_<n>.mat (I, Q, hdr); replay_ddc_udp.m
                                                            %    sends these files with the PRI timing
cfg.io.matFormat = '-v6';                                   % '-v6' fast/uncompressed, '-v7' compressed

cfg.udp.enable       = false;        % stream the DDC output of every PRI over Ethernet/UDP
cfg.udp.host         = '127.0.0.1';
cfg.udp.port         = 50000;
cfg.udp.payloadBytes = 1400;         % DDC data bytes per datagram (int16 I/Q; + 64-byte header fits a 1500-byte
                                     % Ethernet MTU; use up to 65400 with jumbo frames / localhost)

%% Display
cfg.display.enable            = true;
cfg.display.visible           = true;     % false = draw off-screen (figures are still saved)
cfg.display.updateEvery       = 50;       % PPI refresh period [pulses]
cfg.display.maxRange_m        = [];       % [] = cfg.rx.maxRange_m
cfg.display.showTruth         = true;     % draw the true target positions
cfg.display.rangeProfilePulse = [];       % pulse number whose matched filter output is plotted
cfg.display.rangeProfileAxis  = 'range';  % 'range', 'delay' or 'index'
end
