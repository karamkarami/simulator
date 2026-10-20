% RUN_ALL_TESTS  Run every test of the radar simulator and print a summary.
%   The six-target tests take about 10 minutes each in Octave and write about
%   8 GB of raw data each (output/test_<mode>/raw_pulses).

addpath(fileparts(mfilename('fullpath')));
names  = {};
passed = [];

for mode = {'fixed', 'staggered'}
    for pulse = {'short', 'long'}
        names{end+1}  = sprintf('SNR budget (%s PRF, %s pulse)', mode{1}, pulse{1}); %#ok<SAGROW>
        passed(end+1) = test_snr_budget(mode{1}, pulse{1});                      %#ok<SAGROW>
    end
end

names{end+1}  = 'UDP streaming';
passed(end+1) = test_udp_stream();

for mode = {'fixed', 'staggered'}
    names{end+1}  = sprintf('Six targets with clutter (%s PRF)', mode{1});  %#ok<SAGROW>
    passed(end+1) = test_six_targets(mode{1}, true);                       %#ok<SAGROW>
end

fprintf('\n===== Test summary =====\n');
for i = 1:numel(names)
    if passed(i), s = 'PASS'; else, s = 'FAIL'; end
    fprintf('%-45s %s\n', names{i}, s);
end
