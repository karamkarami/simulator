% RUN_ALL_TESTS  Run every test of the radar simulator and print a summary.
%   The six-target tests take about 30 minutes each in Octave. Each writes
%   about 18 GB of raw data (output/test_<mode>/raw_pulses), checks it and
%   deletes it again (call test_six_targets(mode, true) to keep the files).
%   When the C receiver is built (radar_receiver/build/radar_receiver, Linux
%   or macOS), test_c_receiver compares it with the MATLAB processing (about
%   30 minutes per PRF mode).

addpath(fileparts(mfilename('fullpath')));
names  = {};
passed = [];

for mode = {'fixed', 'staggered'}
    for pulse = {'short', 'long'}
        names{end+1}  = sprintf('SNR budget (%s PRF, %s pulse)', mode{1}, pulse{1}); %#ok<SAGROW>
        passed(end+1) = test_snr_budget(mode{1}, pulse{1});                      %#ok<SAGROW>
    end
end

names{end+1}  = 'SNR budget (fixed PRF, long pulse, P4 phase code)';
passed(end+1) = test_snr_budget('fixed', 'long', 'phasecode');
names{end+1}  = 'SNR budget (single element, no beamforming)';
passed(end+1) = test_snr_budget('fixed', 'short', 'lfm', 1);

names{end+1}  = 'UDP streaming of the DDC output';
passed(end+1) = test_udp_stream();

names{end+1}  = 'Interactive PPI (data tips, zoom, pan, filters)';
passed(end+1) = test_ppi_interactive();

names{end+1}  = 'No ghost of the short-pulse echo (same IF)';
passed(end+1) = test_short_pulse_ghost();

for mode = {'fixed', 'staggered'}
    names{end+1}  = sprintf('Six targets with clutter (%s PRF)', mode{1});  %#ok<SAGROW>
    passed(end+1) = test_six_targets(mode{1}, true, false);                %#ok<SAGROW>
end

receiverBin = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'radar_receiver', 'build', 'radar_receiver');
if isunix && exist(receiverBin, 'file')
    for mode = {'fixed', 'staggered'}
        names{end+1}  = sprintf('C receiver vs MATLAB (%s PRF)', mode{1});  %#ok<SAGROW>
        passed(end+1) = test_c_receiver(mode{1}, receiverBin);             %#ok<SAGROW>
    end
else
    fprintf('C receiver not built (%s): test_c_receiver skipped\n', receiverBin);
end

fprintf('\n===== Test summary =====\n');
for i = 1:numel(names)
    if passed(i), s = 'PASS'; else, s = 'FAIL'; end
    fprintf('%-55s %s\n', names{i}, s);
end
