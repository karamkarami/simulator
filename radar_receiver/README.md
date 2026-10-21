# Radar receiver in C (CLion / CMake project)

This program does all the receiver processing of the MATLAB simulator in C.
Its input is the signed 16-bit DDC output of every PRI, which the simulator
(or real hardware using the same packet format) sends over Ethernet/UDP. It
writes the plots to the console and a CSV file. Optionally it also shows them
on a Qt PPI and sends the CFAR detections (or the plots) over a serial port
at 921600 baud.

```
 UDP datagrams (int16 I/Q + header)  -- or a recorded file (--replay)
        |  reassembly of each PRI (udp_input.c)
        v  split in time: short-pulse samples | long-pulse samples
 DBF (21 x 48 matrix, off for 1 element)
 Matched filter (valid convolution; FFT fast convolution when it is faster)
 CPI ring buffer (the pulses of one CPI) -> MTI ring (MTI output of each pulse computed once)
 Doppler: fixed PRF Hamming + FFT 128 | staggered: de-interleave, Hamming, FFT 32, NCI;
          division by the noise gain of every bin
 2-D CA-CFAR -> pruning (local peaks, sidelobe blanking) -> reports -> plot extraction
        |                                   |                         |
        |                                   +-> serial (detections)   +-> CSV, console, Qt PPI,
        |                                                                 serial (plots)
```

The algorithms are the same as in the MATLAB simulator (`../src/*.m`, see
`../README.md`). The C receiver reads every number that MATLAB derives
(DBF matrix, matched filters, range axes, CPI length, MTI weights, Doppler
window and noise gains, CFAR threshold factor, extraction settings) from a
configuration file that `export_receiver_config.m` writes. Both programs
therefore always use the same values.

## Build

Requirements: a C11 compiler (GCC, Clang or MSVC), CMake 3.16 or newer,
OpenMP (optional, for multi-threading), and Qt 5.12+ or Qt 6 Widgets
(optional, for the PPI).

**CLion:** *File > Open* the `radar_receiver` folder; CLion loads
`CMakeLists.txt`. Choose the **Release** profile (*Settings > Build,
Execution, Deployment > CMake*, add a Release profile). Without it the code
is not optimised and runs about 10 times slower. Run configurations are
created for `radar_receiver` and the tests. Set the program arguments in the
run configuration, for example `--config ../config/radar_config.bin --ppi`.
If CMake does not find Qt, add `-DCMAKE_PREFIX_PATH=<Qt>/lib/cmake` (for
example `C:/Qt/6.7.2/mingw_64/lib/cmake`) to the CMake options.

**Command line:**

```sh
cmake -S radar_receiver -B radar_receiver/build -DCMAKE_BUILD_TYPE=Release
cmake --build radar_receiver/build -j
cd radar_receiver/build && ctest --output-on-failure
```

CMake options:

| Option | Default | Meaning |
|---|---|---|
| `RADAR_WITH_QT` | ON | build the PPI (the build continues without Qt when it is not found) |
| `RADAR_OPENMP` | ON | OpenMP threads |
| `RADAR_NATIVE` | ON | optimise for the CPU of the build machine (`-march=native`, MSVC `/arch:AVX2`); switch it off when the program must run on other PCs |
| `RADAR_BUILD_TESTS` | ON | unit tests (`ctest`) |

On Windows, the MinGW toolchain bundled with CLion and MSVC both work.
OpenMP and `-march=native` need no extra setup with MinGW. With MSVC,
`/openmp` is used. The serial port is opened as `COM3` and so on.

## Use with the MATLAB simulator

```matlab
cfg = radar_config;                       % the same cfg for both programs
cfg.udp.enable = true;                    % stream the DDC output
export_receiver_config(cfg);              % -> radar_receiver/config/radar_config.bin
```

```sh
radar_receiver/build/radar_receiver --config radar_receiver/config/radar_config.bin --ppi
```

Then run `run_radar_sim(cfg)` in MATLAB. Export the configuration again
whenever `cfg` changes: the receiver checks that the PRI layout (number of
channels and DDC samples) matches and ignores PRIs that do not.
`tests/test_c_receiver.m` does all of this automatically and compares the
plots of both programs.

## Options

| Option | Meaning |
|---|---|
| `--config FILE` | configuration from `export_receiver_config.m` (default `config/radar_config.bin`) |
| `--port N` | UDP port of the DDC stream (default 50000) |
| `--record FILE` | also store every received datagram in FILE |
| `--replay FILE` | read the datagrams from a recorded FILE instead of the network (as fast as possible) |
| `--ppi` | show the PPI (Qt build) |
| `--ppi-snapshot FILE` | save the PPI as a PNG at the end and close it |
| `--serial DEVICE` | send the CFAR output to a serial port (`/dev/ttyUSB0`, `COM3`, ...) |
| `--baud N` | serial baud rate (default 921600) |
| `--serial-content X` | `detections` (default: the pruned CFAR detections of every CPI update) or `plots` |
| `--threads N` | processing threads (default: number of cores - 1) |
| `--update-every N` | run Doppler, CFAR and extraction every N pulses (default from the configuration: 1) |
| `--plots-csv FILE` | plot list (default `plots.csv`; `''` = none) |
| `--detections-csv FILE` | all pruned CFAR detections (default none) |
| `--idle-timeout S` | stop when no PRI arrived for S seconds after the first one (default 0 = never) |
| `--quiet` | do not print every plot |

Ctrl+C (or closing the PPI window) stops the receiver. The open plots are
closed and the processing times are printed.

## Input: UDP packets

Each datagram is a 64-byte little-endian header followed by int16 (I, Q)
pairs, channel after channel (see "Raw data and UDP format" in
`../README.md`). The receiver reassembles every PRI from its packets in any
order. It ignores duplicates and foreign datagrams, drops PRIs with missing
packets (counted as incomplete), and follows the sender when it restarts.
A jump in the pulse numbers (a PRI lost, or a skipped azimuth sector)
restarts the CPI buffers and closes the open plots, as the MATLAB simulator
does between simulated sectors.

At the full rate of the default radar (48 channels, about 1400 DDC samples
per PRI and 1000 PRIs/s), the stream is 340 MB/s (2.7 Gbit/s). It needs
10 Gbit Ethernet, and the socket receive buffer should be large. The
receiver asks for 32 MB; on Linux, allow it with
`sysctl -w net.core.rmem_max=33554432`.

## Outputs

**Plots** (console, `--plots-csv`, PPI, serial): time, azimuth, range,
elevation (NaN with one element), radial velocity (ambiguous, + = away),
SNR, sub-pulse, number of CPI detections (hits), beam, range cell, Doppler
bin (0-based) and pulse number.

**Serial frames** (little-endian; CRC-16/CCITT-FALSE, polynomial 0x1021,
init 0xFFFF, over the bytes after the sync word up to the CRC):

| Frame | Layout |
|---|---|
| Detections (type 1), 16 + 24 N bytes, N <= 200 per frame | `A5 5A`, u8 type = 1, u8 pulse (0 short, 1 long), u32 pulse number, f32 azimuth [deg], u16 N, then N records: u16 beam, u16 range cell, u16 Doppler bin, u16 reserved, f32 range [m], f32 radial velocity [m/s], f32 elevation [deg], f32 SNR [dB]; u16 CRC |
| Plot (type 2), 42 bytes | `A5 5A`, u8 type = 2, u8 pulse, u32 pulse number, f64 time [s], f32 azimuth, range, elevation, radial velocity, SNR, u16 beam, u16 hits, u16 CRC |

A CPI update with more than 200 detections is sent as several frames with
the same header. The serial output never blocks the processing. Frames go
into a 256 kB buffer; when the port cannot keep up, whole frames are dropped
and counted. 921600 baud carries about 92 kB/s, which is about 3800
detections per second. `tools/decode_serial.py` decodes a port or a capture
file (`decode_serial.py /dev/ttyUSB0`, `decode_serial.py capture.bin --csv`).

**PPI:** range rings, azimuth spokes (0 = North, clockwise) and the antenna
line. Short-pulse plots are cyan and long-pulse plots yellow. They are drawn
at ground range and removed when the antenna passes them again in the next
scan.

## Speed

The processing is written for x86 CPUs with AVX2:
* The data is float and stored as separate real and imaginary arrays.
* The DBF, the matched filter, the Doppler FFTs (radix 2, many signals at
  once, in cache-sized blocks), the MTI and the CFAR (running sums) all
  have inner loops that the compiler vectorises.
* The MTI output of every pulse is computed once and kept in a ring.
* The long matched filter runs as FFT fast convolution when that is
  estimated to be faster than the direct convolution with the available
  threads.
* The stages are split over threads with OpenMP, and all work buffers are
  allocated once.

PERFORMANCE_TABLE

Doppler processing and CFAR on every new pulse, as in MATLAB, need more
time than one PRI. To keep up with the real-time stream, use
`--update-every N` (for example 8 or 16). The fast-time processing (DBF,
matched filter, MTI) still runs for every pulse. Doppler, CFAR and plot
extraction then run every N pulses: 0.29 deg of antenna rotation with N = 8.
The number of CPI detections a plot needs is scaled with N, so the plots
stay the same.

## Tests

| Test | What it checks |
|---|---|
| `core` (`tests/test_core.c`) | batched FFT against a direct DFT; plot extractor; a moving target through the whole chain with the exported configuration (one plot at the right range, velocity and beam); direct and fast-convolution matched filters give the same plot |
| `udp_input` | reassembly with shuffled, duplicate, missing and foreign datagrams; timeouts; record and replay; fuzzing |
| `serial` | frame encoding and CRC; real output to a pseudo-terminal (Linux): exact bytes, 200-detection split, non-blocking drops, close |
| `ppi` | the Qt PPI runs off-screen with antenna updates and plots and saves a snapshot |
| `../tests/test_c_receiver.m` (Octave/MATLAB) | the six-target scenario streamed over UDP to this receiver: the C plots must match the MATLAB plots; the recording is replayed to check that the result is repeatable and to measure the time |
