function results = run_radar_sim(cfg)
%RUN_RADAR_SIM  Simulation of a rotating S-band search radar with elevation DBF.
%
%   results = RUN_RADAR_SIM()    runs the simulation with radar_config().
%   results = RUN_RADAR_SIM(cfg) runs it with a modified configuration.
%
%   Processing chain of every PRI (one long + one short sub-pulse):
%     1. Signal generation : target echoes + ground clutter + noise at the IF,
%                            48 elements, fs = 120 MHz, 16-bit ADC (12.7 ENOB)
%     2. Recording         : raw ADC data of every PRI -> <rawFolder>/pulse_<n>.mat
%                            (optional UDP streaming)
%     3. DDC               : two stages, decimation 10 and 6 (valid convolutions)
%     4. DBF               : 48 elements -> 21 elevation beams
%     5. Matched filter    : per beam (valid convolution, Taylor window for LFM)
%     6. CPI buffers       : circular buffers holding the pulses of one beamwidth
%     7. MTI + Doppler     : 3-pulse canceller + 32-point FFT (fixed PRF), or
%                            stagger MTI + de-interleaving + FFT + NCI (staggered PRF)
%     8. 2-D CFAR          : on every range-Doppler map
%     9. Pruning           : local peaks only, then azimuth peak picking (plots)
%    10. PPI display       : updated while the radar is running
%
%   results.plots holds the extracted plots (time, azimuth, range, elevation,
%   ambiguous radial velocity, SNR, pulse type).

if nargin < 1, cfg = radar_config(); end
root = fileparts(mfilename('fullpath'));
addpath(root, fullfile(root, 'src'));
rng(cfg.sim.seed);
tStart = tic;

%% Setup
sys   = setup_radar(cfg);
sched = pulse_schedule(sys);
names = {'short', 'long'};
if cfg.sim.verbose
    print_summary(sys, sched);
    results.budget = print_snr_budget(sys, cfg.targets, 0);
end

clut = [];
if cfg.clutter.enable
    clut = clutter_init(sys);
end
if cfg.io.saveRaw && ~exist(cfg.io.rawFolder, 'dir')
    mkdir(cfg.io.rawFolder);
end
if ~exist(cfg.io.outFolder, 'dir')
    mkdir(cfg.io.outFolder);
end
udp = [];
if cfg.udp.enable
    udp = udp_sender_open(cfg.udp);
end
ppi = [];
if cfg.display.enable
    ppi = ppi_init(sys);
end

for i = 1:2
    p = names{i};
    buf.(p) = cpi_buffer_init(sys.nBeams, sys.wave.(p).nRange, sys.cpi.length);
    ext.(p) = plot_extractor_init(sys.dop.nfft);
end
plots = [];
nProcessed = 0;
wasActive = false;
% Plot extraction options: minimum azimuth extent -> minimum number of CPI detections
extractOpt = cfg.extract;
extractOpt.minHits = max(1, ceil(cfg.extract.minExtent_deg / ...
                     (sys.ant.rate_dps * mean(sys.pri) * cfg.proc.updateEvery)));
extractOpt.flush = false;
flushOpt = extractOpt;
flushOpt.flush = true;

%% Pulse loop
for k = 1:sched.nPulses
    t  = sched.time_s(k);
    az = sched.az_deg(k);

    % Skipped sector: restart the CPI buffers and close the open plots
    if ~sched.active(k)
        if wasActive
            for i = 1:2
                p = names{i};
                buf.(p) = cpi_buffer_init(sys.nBeams, sys.wave.(p).nRange, sys.cpi.length);
                [ext.(p), newPlots] = plot_extractor(ext.(p), [], flushOpt);
                plots = [plots, newPlots]; %#ok<AGROW>
            end
        end
        wasActive = false;
        continue
    end
    wasActive = true;

    % 1) Received signal of all elements and ADC
    clutterIf = [];
    if ~isempty(clut)
        clutterIf = clutter_signal(clut, sys, az, t);
    end
    adcData = adc_sample(generate_rx_signal(sys, cfg.targets, t, az, clutterIf), sys);

    % 2) Raw data recording and streaming
    hdr = pulse_header(sys, sched, k);
    if cfg.io.saveRaw
        save_pulse_file(cfg.io.rawFolder, cfg.io.matFormat, hdr, adcData);
    end
    if ~isempty(udp)
        udp_send_pulse(udp, hdr, adcData);
    end

    % 3-9) Processing chain, separately for the short and the long sub-pulse
    for i = 1:2
        p  = names{i};
        wv = sys.wave.(p);
        bb = ddc(adcData, sys, wv);                    % [nSamples x nElements]
        bm = digital_beamforming(bb, sys.dbf.W);        % [nBeams x nDdc]
        rp = matched_filter(bm, wv);                    % [nBeams x nRange]
        if k == cfg.display.rangeProfilePulse
            plot_range_profile(rp, wv, cfg.display.rangeProfileAxis);
        end

        buf.(p) = cpi_buffer_push(buf.(p), rp, sched.staggerIdx(k), az, t);
        if buf.(p).count == sys.cpi.length && mod(buf.(p).nPushed, cfg.proc.updateEvery) == 0
            [cpi, info] = cpi_buffer_read(buf.(p));
            y  = mti_filter(cpi, sys.mtiCoeffs, info.staggerIdx);
            rd = doppler_process(y, sys);               % [nBeams x nRange x nfft]
            [isDet, noise] = cfar_2d(rd, sys);
            dets = prune_peaks(isDet, rd, noise, sys, wv);
            reports = detection_reports(dets, rd, sys, wv, info);
            [ext.(p), newPlots] = plot_extractor(ext.(p), reports, extractOpt);
            plots = [plots, newPlots]; %#ok<AGROW>
        end
    end
    nProcessed = nProcessed + 1;

    % 10) Display
    if ~isempty(ppi) && mod(nProcessed, cfg.display.updateEvery) == 0
        ppi_update(ppi, sys, az, t, sched.scan(k), plots);
    end
    if cfg.sim.verbose && mod(nProcessed, 500) == 0
        fprintf('  pulse %6d (az %6.1f deg): %d plots so far, %.0f s elapsed\n', k, az, numel(plots), toc(tStart));
    end
end

%% Close the remaining plots
for i = 1:2
    p = names{i};
    [ext.(p), newPlots] = plot_extractor(ext.(p), [], flushOpt);
    plots = [plots, newPlots]; %#ok<AGROW>
end
if ~isempty(ppi)
    ppi_update(ppi, sys, sched.az_deg(end), sched.time_s(end), sched.scan(end), plots);
    print(ppi.fig, fullfile(cfg.io.outFolder, 'ppi.png'), '-dpng');
end
if ~isempty(udp)
    udp_sender_close(udp);
end

results.plots      = plots;
results.sys        = sys;
results.nProcessed = nProcessed;
results.runTime_s  = toc(tStart);
if cfg.sim.verbose
    print_plots(plots);
    fprintf('%d PRIs simulated and processed in %.1f s (%.1f ms per PRI)\n', nProcessed, ...
        results.runTime_s, 1e3*results.runTime_s/max(nProcessed, 1));
end
end

function hdr = pulse_header(sys, sched, k)
% Information stored with the raw data of PRI k.
hdr.pulseNumber    = k;
hdr.scanNumber     = sched.scan(k);
hdr.time_s         = sched.time_s(k);
hdr.antennaAz_deg  = sched.az_deg(k);
hdr.pri_s          = sched.pri_s(k);
hdr.staggerIndex   = sched.staggerIdx(k);
hdr.fs_Hz          = sys.fs;
hdr.rxStart_s      = sys.rxStart_s;
hdr.fc_Hz          = sys.cfg.fc;
hdr.longTx_s       = [sys.wave.long.tTx_s,  sys.wave.long.width_s];
hdr.shortTx_s      = [sys.wave.short.tTx_s, sys.wave.short.width_s];
hdr.longIF_Hz      = sys.wave.long.fIF_Hz;
hdr.shortIF_Hz     = sys.wave.short.fIF_Hz;
end

function print_summary(sys, sched)
% Short description of the simulated radar.
cfg = sys.cfg;
fprintf('\n=== Rotating search radar simulation ===\n');
fprintf('fc = %.2f GHz, Pt = %.0f kW, NF = %.1f dB, %d elements -> %d beams\n', cfg.fc/1e9, ...
    cfg.tx.peakPower_W/1e3, cfg.rx.noiseFigure_dB, sys.nCh, sys.nBeams);
fprintf('PRF mode: %s, mean PRI = %.0f us, %d rpm, azimuth beamwidth %.2f deg\n', cfg.prf.mode, ...
    mean(sys.pri)*1e6, cfg.ant.rpm, sys.ant.azBeamwidth_deg);
fprintf('Receive window: %.1f ... %.1f us (%d ADC samples, %.1f MB per PRI)\n', sys.rxStart_s*1e6, ...
    sys.tAdc(end)*1e6, sys.nAdc, sys.nAdc*sys.nCh*2/2^20);
fprintf('Long  pulse: %s %3.0f us, IF %.0f MHz, range %.2f ... %.2f km (%d cells)\n', sys.wave.long.type, ...
    sys.wave.long.width_s*1e6, sys.wave.long.fIF_Hz/1e6, sys.wave.long.range_m([1 end])/1e3, sys.wave.long.nRange);
fprintf('Short pulse: %s %3.0f us, IF %.0f MHz, range %.2f ... %.2f km (%d cells)\n', sys.wave.short.type, ...
    sys.wave.short.width_s*1e6, sys.wave.short.fIF_Hz/1e6, sys.wave.short.range_m([1 end])/1e3, sys.wave.short.nRange);
fprintf('CPI: %d pulses = %d row(s) x %d samples, FFT %d, Doppler bin %.2f m/s, unambiguous +-%.1f m/s\n', ...
    sys.cpi.length, sys.cpi.nRows, sys.cpi.rowLen, sys.dop.nfft, abs(diff(sys.dop.v_mps(1:2))), ...
    sys.lambda/(4*sys.dop.rowPri_s));
fprintf('CFAR: Pfa = %.0e, %d training cells, %.1f looks, alpha = %.2f (%.1f dB)\n', cfg.cfar.pfa, ...
    sys.cfar.nTrain, sys.dop.nLooks, sys.cfar.alpha, 10*log10(sys.cfar.alpha));
fprintf('Simulated pulses: %d of %d (%.1f s of radar time)\n', sum(sched.active), sched.nPulses, sched.time_s(end));
end

function print_plots(plots)
% List the extracted plots.
fprintf('\nExtracted plots: %d\n', numel(plots));
fprintf('%8s %8s %9s %7s %8s %7s %6s %5s\n', 'time[s]', 'az[deg]', 'R[km]', 'el[deg]', 'vr[m/s]', ...
    'SNR[dB]', 'pulse', 'hits');
for i = 1:numel(plots)
    p = plots(i);
    fprintf('%8.3f %8.2f %9.3f %7.2f %8.2f %7.1f %6s %5d\n', p.time_s, p.az_deg, p.range_m/1e3, ...
        p.el_deg, p.vr_mps, p.snr_dB, p.pulse, p.hits);
end
end
