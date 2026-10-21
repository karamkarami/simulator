# Rotating search radar simulator (MATLAB)

A simple, modular MATLAB simulation of a 3 GHz rotating search radar. The
receive antenna is a column of 48 elements (selectable) with digital
beamforming (DBF) in elevation, and every PRI holds a short and a long
sub-pulse. Everything from the 120 MHz ADC samples on the 30 MHz IF down to
the PPI display is simulated. Each processing block is a short function in
`src/`. All parameters are in one file, `radar_config.m`.

The code needs no MATLAB toolbox (UDP streaming uses `udpport` or Java). It
was tested with GNU Octave 8.4. The folder `radar_receiver/` holds a C
version of the whole receiver processing (CMake project for CLion). It
processes the DDC output that the simulator streams over UDP and has an
optional Qt PPI and serial output (see `radar_receiver/README.md`).

```
 PRI:  short pulse 25 us at 0 | receive (near range) | long pulse 150 us at 200 us | receive (far range)
 per target: position from its motion -> power from the radar equation (Swerling 1)
 echoes + clutter + Tx leakage (40 dB isolation) + thermal noise, 30 MHz IF, 48 elements
 x receiver gain (60 dB) -> ADC 120 MHz, 16 bit, ENOB 12.7 -> int16
        |                       \--> raw data: output/raw_pulses/pulse_<n>.mat
        v
 DDC: x*cos(30 MHz) | x*(-sin(30 MHz)) -> FIR h1, /10 -> FIR h2, /6
      (int16 coefficients with 13 fractional bits, time-domain sliding valid convolutions,
       every stage rounded and saturated to int16)                     -> int16 [nDdc x 48]
        |                       \--> optional UDP: int16 DDC output + pulse/scan/azimuth header
        |                             (-> radar_receiver/, the same processing in C)
        v  split in time: short-pulse samples | long-pulse samples
 DBF (21 x 48 matrix, off for 1 element)                   -> [21 beams x nSamples] per sub-pulse
 Matched filter (valid convolution, Taylor for LFM)        -> [21 beams x nRange]   (range/delay axes)
 CPI buffer (circular, the pulses of one CPI)              -> [21 x nRange x 114], short and long
 MTI  fixed PRF: 3-pulse canceller | staggered: weights from file
 Doppler fixed: Hamming + FFT 128 (next power of 2 of 114) | staggered: de-interleave, Hamming, FFT 32, NCI
 2-D CA-CFAR on every range-Doppler map
 Pruning (local peaks + sidelobe blanking) -> plot extraction (azimuth peak picking)
 PPI display (updated while the antenna turns)
```

## Quick start

```matlab
results = run_radar_sim;            % default configuration (full 360 deg scan)

cfg = radar_config;                 % or change parameters first
cfg.prf.mode = 'staggered';
cfg.clutter.enable = false;
cfg.ant.nElements = 1;              % a single antenna: no beamforming
cfg.sim.sectors_deg = [20 45];      % only simulate this azimuth sector
results = run_radar_sim(cfg);

cd tests
run_all_tests                       % SNR budget, UDP, ghost and six-target tests
```

C receiver (see `radar_receiver/README.md`):

```matlab
cfg = radar_config;  cfg.udp.enable = true;
export_receiver_config(cfg);        % configuration file for the C receiver
% start radar_receiver/build/radar_receiver --ppi, then:
results = run_radar_sim(cfg);
```

`make_input_files.m` creates the default input files in `inputs/`. The
simulator calls it automatically when a file is missing, or when a default
file was built for other values (for example after you change
`cfg.ant.nElements` or `cfg.ant.tilt_deg`). Replace a file with your own data
(same variable names) to simulate another antenna, filter set or stagger.

## Files

| File | Purpose |
|---|---|
| `radar_config.m` | all user parameters |
| `export_receiver_config.m` | writes the configuration of the C receiver (`radar_receiver/config/radar_config.bin`) |
| `run_radar_sim.m` | main loop: generation, recording, DDC/streaming, processing, display |
| `make_input_files.m` | default input .mat files |
| `src/setup_radar.m` | derives the PRI timeline, NCO, DDC segments, matched filters, axes, CPI length, FFT length, CFAR factor ... |
| `src/generate_rx_signal.m`, `clutter_init.m`, `clutter_signal.m`, `adc_sample.m` | signal generation |
| `src/target_rcs.m`, `draw_target_rcs.m` | Swerling 1 target RCS (a new value every scan) |
| `src/save_pulse_file.m` | raw ADC data of a PRI to a .mat file |
| `src/ddc.m`, `fir_decimate_valid.m` | two-stage DDC with I/Q mixer |
| `src/udp_sender_open.m`, `udp_send_pulse.m`, `udp_sender_close.m`, `udp_receiver_open.m`, `udp_receive_pulse.m` | UDP streaming of the DDC output (sender and reference receiver) |
| `src/digital_beamforming.m` | DBF |
| `src/matched_filter.m`, `plot_range_profile.m` | pulse compression, range/delay/index plots |
| `src/cpi_buffer_*.m` | circular CPI buffers |
| `src/mti_filter.m`, `doppler_process.m` | clutter filter, Doppler processing |
| `src/cfar_2d.m`, `cfar_alpha.m` | 2-D CA-CFAR |
| `src/prune_peaks.m`, `detection_reports.m`, `plot_extractor*.m` | pruning and plot extraction |
| `src/ppi_init.m`, `ppi_update.m` | PPI display |
| `src/snr_budget.m`, `print_snr_budget.m` | expected SNR after every block |
| `src/six_targets_scenario.m`, `make_target.m`, `target_geometry.m`, `illumination_time.m`, `target_sectors.m` | targets |
| `tests/` | tests (see below) |
| `radar_receiver/` | the receiver processing in C (CMake/CLion, Qt PPI, serial output) |

## Input files (`inputs/`)

| File | Variables |
|---|---|
| `antenna_pattern.mat` | `elementSpacing_m` (0.66 wavelength), `az_deg`, `azGain_dB` (one-way, 4.1 deg), `elRel_deg`, `elElemGain_dB` (element pattern vs. angle from broadside, 70 deg), `txElGain_dB`, `rxChannelGain_dBi`, `txGain_dBi` |
| `ddc_filters.mat` | `h1` (stage 1 at 120 MHz, 70 taps), `h2` (stage 2 at 12 MHz, 66 taps): int16 with 13 fractional bits (coefficient = value / 8192), DC gain about 1 |
| `dbf_matrix.mat` | `W` [21 x nElements] (beam = W * x), `beamEl_deg` |
| `stagger_prf.mat` | `priVector_s` [K x 1], `mtiCoeffs` [K x L] (row k = weights of an MTI output whose first pulse has stagger position k; one row = fixed weights) |

The default files also hold `madeFor`, the configuration values they were
built for. User files simply leave it out.

## PRI timeline and range coverage

| Time in the PRI | Event |
|---|---|
| 0 - 25 us | short pulse transmitted (`cfg.wave.short.tTx_s`); Tx leakage in the receiver |
| 25 - 200 us | short-pulse receive interval |
| 200 - 350 us | long pulse transmitted (`cfg.wave.long.tTx_s`); Tx leakage in the receiver |
| 350 - 896 us | long-pulse receive interval, up to `cfg.rx.maxRange_m` (80 km) |

The transmit times are parameters. Both sub-pulses use LFM (1 MHz) or a P4
phase code (1 us chips) on the same 30 MHz IF and are separated only in
time. The ADC records one continuous block from the start of the PRI to the
end of the long window, Tx leakage included. The DDC processes the whole
block, and its output is then split: a DDC sample belongs to a sub-pulse
when its whole filter input lies inside that sub-pulse's receive interval.
The leakage samples are therefore never processed. The short pulse covers
4.2 - 24.5 km and the long pulse 22.9 - 81.4 km (the long pulse is blind up
to c x 150 us / 2 = 22.5 km). The CFAR does not test the outermost 11 range
cells of each window, so the short pulse hands over to the long pulse at
23.7 km. The latest long-pulse start, and with it the instrumented range, is
limited by the shortest PRI; `setup_radar` stops with an error when a
receive interval is too short for the filters and the CFAR window.

Because both sub-pulses share the IF, the short-pulse echo of a target
beyond 52.5 km also falls into the long-pulse window. The long-pulse matched
filter does not match it, so it spreads over about 10 km of range, 30 - 52 km
in front of the target, at up to about 24 dB below the target. The 2-D CFAR
adapts to it like a local rise of the noise floor, so it causes no false
plots (see `test_short_pulse_ghost`). Its only effect is lower sensitivity
there, at the same beam and Doppler, behind very strong far targets.

## Model summary

* **Geometry:** antenna 6 m above flat ground, tilted 11 deg up, turning
  clockwise at 6 rpm. Targets fly with a constant velocity (`make_target.m`).
  Azimuth is measured clockwise from North.
* **Received signal** (`generate_rx_signal.m`, `adc_sample.m`), built in
  this order for every PRI:
  1. For every target, its position from its motion and the received power
     from the radar equation (20 kW, antenna patterns, 3 dB losses, RCS of
     the current scan). The echo of both sub-pulses on the 30 MHz IF has the
     element phases of the tilted column and the Doppler inside each pulse
     (stop-and-hop between pulses).
  2. Clutter.
  3. Transmitter leakage while a pulse is transmitted: Pt - 40 dB
     isolation, in every channel.
  4. Thermal noise, k T0 F over fs/2 (NF 4 dB).
  5. The receiver gain of 60 dB.
  6. The ADC: full scale +10 dBm (a 2 Vpp sine into 50 ohm), white ADC noise
     for an ENOB of 12.7 added at the quantizer input, then rounding and
     clipping to int16.

  Signals are in sqrt(W) at the receiver input. Thermal noise reaches the
  ADC at 180 LSB rms (-42 dBFS), 36 dB above the ADC noise (2.8 LSB rms).
  The leakage is +93 dBm at the ADC and saturates it while a pulse is
  transmitted. It never reaches the processing (see the PRI timeline). Very
  strong near clutter can clip the ADC as well, because there is no STC.
* **Targets:** Swerling 1. Every target gets a new exponentially
  distributed RCS (mean `rcs_m2`) every scan (`cfg.sim.swerling = 0` for
  constant RCS).
* **Clutter:** constant sigma0 (-30 dB) Rayleigh cells (0.25 deg x 15 m) from
  1 km to the radar horizon, with internal motion (0.1 m/s). Scan modulation
  comes from the two-way azimuth pattern. The echo is convolved with both
  transmitted IF pulses, so partly received (eclipsed) echoes are included.
* **DDC (fixed point):** the real line is `x .* cos(2*pi*30MHz*t)` and the
  imaginary line `x .* (-sin(2*pi*30MHz*t))` (NCO values with 15
  fractional bits). Each line then goes through h1 with decimation by 10
  and h2 with decimation by 6. Both are valid convolutions in the time
  domain: the filter slides over the signal and is evaluated only at the
  kept output samples (`fir_decimate_valid.m`, equal to
  `conv(x,h,'valid')` followed by `y(1:D:end)`; no FFT). The filter
  coefficients are integers with 13 fractional bits. Each stage shifts its
  accumulator right by 13 bits (`cfg.ddc.shift`, unit gain), rounds and
  saturates it to int16. The DDC output is therefore signed 16-bit, with
  the thermal noise at about 16 LSB rms per I/Q component.
* **Elements and DBF:** the antenna has a fixed 4.1 deg beam in azimuth
  and a column of `cfg.ant.nElements` elements (default 48) spaced 0.66
  wavelength in elevation. The default DBF matrix forms 21 beams equally
  spaced in sine space with a Taylor taper. With 1 element the DBF block is
  skipped, the single channel goes straight to the matched filter, and
  plots carry no elevation (NaN).
* **Axes:** the delay and range of every matched-filter output cell are
  computed in `setup_radar.m` (`wv.delay_s`, `wv.range_m`). They follow from
  the ADC sample times, the valid convolutions of the DDC and matched
  filters, and the transmit time of the sub-pulse.
  `plot_range_profile(rp, wv, 'range'|'delay'|'index')` plots against any
  of them. Set `cfg.display.rangeProfilePulse = <pulse number>` to plot it
  during a run.
* **CPI and buffers:** the CPI is the time a target stays inside the 3 dB
  azimuth beam: 4.1 deg / 36 deg/s = 113.9 ms. Each buffer holds the pulses
  that arrive in that time, `round(CPI / PRI)`: 114 with the fixed PRF. With
  a staggered PRF the number is reduced to a multiple of the number of PRIs
  in the stagger: 116 with the default 4-PRI stagger (mean PRI 980 us). The
  buffers are circular, and MTI, Doppler and CFAR run on every new pulse
  (`cfg.proc.updateEvery`).
* **Doppler:** fixed PRF: 3-pulse canceller (112 outputs), Hamming window
  and one FFT whose length is the next power of 2 of the number of pulses
  (114 -> 128). Staggered PRF: the stagger MTI weights from the file, then
  de-interleaving into 4 rows of 28 samples. Each row is sampled every
  sum(PRI); it gets a Hamming window and a 32-point FFT (the next power of 2
  of 28), and |X|^2 is averaged over the rows.
* **CFAR:** 2-D CA-CFAR, 8 + 3 (guard) range cells and 4 + 3 Doppler cells on
  each side, Pfa = 1e-8. The MTI colours the noise ([1 -2 1] gives +4 dB
  near PRF/2), so every Doppler bin is first divided by its noise gain and
  the CFAR sees a flat noise floor. The threshold factor uses the effective
  number of independent rows: 1 with a fixed PRF, 1.5 with the default
  stagger, because the de-interleaved rows are built from overlapping MTI
  outputs. The bins within one Doppler window mainlobe of zero (+-2 bins: the
  MTI notch, where the clutter residue is) are not tested.
* **Pruning and plots:** only local maxima in beam x range x Doppler are
  kept. Range, elevation and Doppler sidelobes more than 20 dB below a
  stronger peak are removed. Detections of the same target in successive
  (sliding) CPIs are grouped, and the strongest one becomes a plot. A plot
  needs detections over at least 2 deg of antenna rotation, half a beamwidth
  (`cfg.extract.minExtent_deg`). Noise false alarms last about 1 deg, while
  a target is detected over about 6-8 deg. Elevation is interpolated between adjacent
  beams. The plot SNR is measured against the noise floor of the beam. The
  radial velocity is ambiguous (+-25 m/s with the fixed PRF, +-6.4 m/s with
  the default stagger).

## Raw data and UDP format

Every simulated PRI is saved as `output/raw_pulses/pulse_<number>.mat`
(`cfg.io.rawFolder`). `adcData` is int16 [nSamples x nElements]. `hdr`
holds the pulse number, scan number, time, antenna azimuth, PRI, stagger
index, fs, IF, carrier, record start, pulse timing, DDC rate and the DDC
samples of each sub-pulse. At 80 km instrumented range one PRI is about
9.5 MB, so a full scan (10 000 PRIs) is about 95 GB. Use
`cfg.sim.sectors_deg` or a smaller `cfg.rx.maxRange_m` for long runs.

With `cfg.udp.enable = true`, the **DDC output** of every PRI (signed
16-bit I/Q, [nDdc x nElements], both sub-pulses) is sent to
`cfg.udp.host:port`. Each datagram has a 64-byte little-endian header
followed by int16 (I, Q) pairs, channel after channel. The default payload
of 1400 bytes keeps each datagram within a standard 1500-byte Ethernet MTU.
One PRI of the default radar is 48 x 1781 x 4 bytes = 342 kB in 245
datagrams.

| Bytes | Type | Field |
|---|---|---|
| 0-3 | uint32 | magic 0x52444443 ('RDDC') |
| 4-7, 8-11 | uint32 | pulse number, scan number |
| 12-15 | single | antenna azimuth [deg] |
| 16-17, 18-19 | uint16 | packet index (0-based), packet count |
| 20-21, 22-23 | uint16 | channels, stagger index |
| 24-27, 28-31 | uint32 | samples per channel, first sample in this packet (0-based) |
| 32-39 | double | PRI start time [s] |
| 40-43, 44-47 | single | DDC sample rate [Hz], time of the first DDC sample in the PRI [s] |
| 48-63 | uint32 | short pulse first/last DDC sample, long pulse first/last DDC sample (1-based) |

`src/udp_receive_pulse.m` is a reference receiver that reassembles one PRI.

## SNR budget

`print_snr_budget` (called at the start of every run) prints the expected
SNR of every target after the ADC, DDC, DBF, matched filter, MTI and Doppler
blocks. Each target is evaluated at the moment the beam crosses it. The
table also gives the SNR needed for Pd = 0.9 (Albersheim) and the margin.
Where a closed form would be approximate, the budget runs the real
processing functions on a noise-free target:
* the matched-filter loss: the pulse, with its Doppler shift, goes through
  the DDC filters and the matched filter. This covers the Taylor window,
  the DDC filtering of phase-code chips, and the Doppler mismatch.
* the slow-time gain: the MTI and Doppler functions run on the target and
  on white noise. This covers the MTI response, the Hamming loss, Doppler
  straddle and the azimuth beam shape over the CPI.

`tests/test_snr_budget.m` measures the SNR after every block of the
simulated chain (signal-only and noise-only runs) and compares it with the
budget.

## Tests

| Test | What it checks |
|---|---|
| `test_snr_budget(mode, pulse, type, nElements)` | measured vs. expected SNR after each block (LFM or P4, 1 or 48 elements) |
| `test_udp_stream` | the DDC output of a PRI sent over UDP to localhost is received bit-exact, with its header |
| `test_short_pulse_ghost` | a strong target at 65 km is detected, and its short-pulse echo in the long window gives no ghost plot |
| `test_six_targets(mode, clutterOn, keepRaw)` | 6 targets (0.5 - 10 m^2 mean RCS, Swerling 1, 8 - 48 km, up to 290 m/s, around the radar) are detected with the right range, azimuth, elevation and ambiguous velocity, with at most 2 false plots; one raw file per PRI, named by the pulse number |
| `test_c_receiver(mode)` | the six-target scenario streamed over UDP to the C receiver (`radar_receiver/`): its plots must match the MATLAB plots; the recorded stream is replayed to check repeatability and speed |

`tests/run_all_tests.m` runs all of them for the fixed and the staggered
PRF. The six-target test simulates only the azimuth sectors around the
targets, about 2400 PRIs. Its raw files are about 18 GB per PRF mode;
`run_all_tests` deletes them after checking them.

## Assumptions and limits

* Flat earth for target geometry. The clutter extends to the 4/3-earth
  radar horizon, and all clutter arrives from one elevation (0 deg).
* The transmitter uses the wide element pattern in elevation. Every
  receive element sees the full azimuth pattern.
* Swerling 1 targets (constant RCS within a scan), no multipath, no echoes
  from the previous PRI.
* LFM range-Doppler coupling is not corrected (about 130 m at 290 m/s with the
  150 us pulse), because the measured velocity is ambiguous.
* Signal generation uses single precision. That is enough for a 16-bit ADC
  and about twice as fast. The DDC computes with integers (exact in
  double). The processing after the DDC uses single precision, as the C
  receiver does.
