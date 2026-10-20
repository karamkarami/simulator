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
%   Quantities that come from input .mat files (antenna patterns, number of
%   receive elements, DDC filters, DBF matrix, staggered PRIs and MTI weights)
%   are listed under "Input files". Default files are created by
%   make_input_files.m when they are missing.

root = fileparts(mfilename('fullpath'));
addpath(fullfile(root, 'src'));

%% Physical constants
cfg.const.c  = 299792458;       % speed of light [m/s]
cfg.const.kB = 1.380649e-23;    % Boltzmann constant [J/K]
cfg.const.T0 = 290;             % reference noise temperature [K]

%% Input files
cfg.files.antenna = fullfile(root, 'inputs', 'antenna_pattern.mat'); % patterns, nElements, spacing, gains
cfg.files.ddc     = fullfile(root, 'inputs', 'ddc_filters.mat');     % low-pass filters of the 2 DDC stages
cfg.files.dbf     = fullfile(root, 'inputs', 'dbf_matrix.mat');      % elevation beamforming matrix
cfg.files.stagger = fullfile(root, 'inputs', 'stagger_prf.mat');     % staggered PRIs + MTI weights

%% Radar site and antenna
cfg.fc           = 3e9;   % center frequency [Hz]
cfg.ant.height_m = 6;     % antenna height above the ground [m]
cfg.ant.tilt_deg = 11;    % array tilt above the horizon [deg]
cfg.ant.rpm      = 6;     % rotation speed [rev/min]
cfg.ant.az0_deg  = 0;     % antenna azimuth at t = 0, clockwise from North [deg]

%% Transmitter and receiver
cfg.tx.peakPower_W     = 20e3;   % peak transmit power [W]
cfg.rx.noiseFigure_dB  = 4;      % receiver noise figure [dB]
cfg.rx.losses_dB       = 3;      % system losses not modelled elsewhere (RF, radome, ...) [dB]
cfg.rx.maxRange_m      = 100e3;  % instrumented range: sets the length of the receive window [m]
cfg.rx.addThermalNoise = true;   % false only for debugging / SNR measurements

%% ADC
cfg.adc.fs_Hz            = 120e6;  % sampling frequency [Hz]
cfg.adc.nBits            = 16;     % ADC word length [bit]
cfg.adc.enob             = 12.7;   % effective number of bits
cfg.adc.noiseAboveAdc_dB = 10;     % receiver gain: thermal noise level above the ADC noise floor [dB]
cfg.adc.addAdcNoise      = true;   % false only for debugging / SNR measurements
cfg.adc.fIFcenter_Hz     = 30e6;   % IF of the center frequency fc (receiver LO = fc - 30 MHz)

%% Waveforms: every PRI holds a long and a short sub-pulse
% The long pulse is transmitted first and the short pulse right after it.
% Both use a different IF (RF = fc + fIF - fIFcenter) so the two DDC channels separate them.
cfg.wave.long.type         = 'lfm';     % 'lfm' or 'phasecode' (P4 polyphase code)
cfg.wave.long.width_s      = 150e-6;    % pulse width [s]
cfg.wave.long.bandwidth_Hz = 1e6;       % LFM sweep bandwidth [Hz]
cfg.wave.long.chip_s       = 1e-6;      % phase-code chip length [s]
cfg.wave.long.fIF_Hz       = 22e6;      % IF carrier at the ADC [Hz]
cfg.wave.long.procRange_m  = [0 inf];   % range interval kept after the matched filter [m]

cfg.wave.short.type         = 'lfm';
cfg.wave.short.width_s      = 25e-6;
cfg.wave.short.bandwidth_Hz = 1e6;
cfg.wave.short.chip_s       = 1e-6;
cfg.wave.short.fIF_Hz       = 38e6;
cfg.wave.short.procRange_m  = [0 30e3]; % short pulse only fills the blind zone of the long pulse

%% PRF
cfg.prf.mode  = 'fixed';   % 'fixed' or 'staggered' (staggered PRIs come from cfg.files.stagger)
cfg.prf.pri_s = 1e-3;      % PRI in fixed mode [s]

%% DDC (filters come from cfg.files.ddc)
cfg.ddc.decimation = [10 6];   % decimation factor of stage 1 and stage 2 (120 -> 12 -> 2 MHz)

%% Matched filter and Doppler processing
cfg.proc.taylorNbar   = 4;     % Taylor window (LFM only): number of nearly constant sidelobes
cfg.proc.taylorSll_dB = -35;   % Taylor window peak sidelobe level [dB]
cfg.proc.cpiLength    = [];    % pulses per CPI buffer ([] = pulses inside the 3 dB azimuth beam)
cfg.proc.nfft         = 32;    % Doppler FFT length
cfg.proc.updateEvery  = 1;     % run MTI/Doppler/CFAR every N new pulses (1 = every pulse)

%% 2-D CA-CFAR on every range-Doppler map
cfg.cfar.pfa          = 1e-7;  % design probability of false alarm
cfg.cfar.guardRange   = 3;     % guard cells on each side in range
cfg.cfar.trainRange   = 8;     % training cells on each side in range
cfg.cfar.guardDoppler = 2;     % guard cells on each side in Doppler
cfg.cfar.trainDoppler = 2;     % training cells on each side in Doppler
cfg.cfar.notchBins    = 1;     % Doppler bins on each side of zero that are not tested (MTI notch)

%% Pruning of the CFAR detections
cfg.prune.sidelobeMargin_dB = 20;  % a peak this much weaker than a peak on the same range, beam
                                   % or Doppler line is treated as its sidelobe and removed

%% Plot extraction (azimuth peak picking over successive CPIs)
cfg.extract.minExtent_deg = 1; % minimum azimuth extent of the detections of a plot [deg]
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
cfg.sim.verbose     = true;

%% Targets (see make_target.m)
cfg.targets = six_targets_scenario();

%% Recording and streaming of the raw ADC data
cfg.io.outFolder = fullfile(root, 'output');
cfg.io.saveRaw   = true;                                    % one .mat file per PRI
cfg.io.rawFolder = fullfile(root, 'output', 'raw_pulses');
cfg.io.matFormat = '-v6';                                   % '-v6' fast/uncompressed, '-v7' compressed

cfg.udp.enable       = false;        % stream the raw data of every PRI over Ethernet/UDP
cfg.udp.host         = '127.0.0.1';
cfg.udp.port         = 50000;
cfg.udp.payloadBytes = 8192;         % ADC bytes per datagram (a 40-byte header is added)

%% Display
cfg.display.enable            = true;
cfg.display.visible           = true;     % false = draw off-screen (figures are still saved)
cfg.display.updateEvery       = 50;       % PPI refresh period [pulses]
cfg.display.maxRange_m        = [];       % [] = cfg.rx.maxRange_m
cfg.display.showTruth         = true;     % draw the true target positions
cfg.display.rangeProfilePulse = [];       % pulse number whose matched filter output is plotted
cfg.display.rangeProfileAxis  = 'range';  % 'range', 'delay' or 'index'
end
