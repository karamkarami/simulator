function export_receiver_config(cfg, fileName)
%EXPORT_RECEIVER_CONFIG  Write the receiver configuration for the C receiver (radar_receiver/).
%
%   EXPORT_RECEIVER_CONFIG(cfg, fileName) runs setup_radar(cfg) and writes
%   everything the C receiver needs to process the DDC output that the
%   simulator streams over UDP (DBF matrix, matched filters, range axes, CPI
%   length, MTI weights, Doppler window and noise gains, CFAR factor, plot
%   extraction settings). Export again whenever cfg changes, before starting
%   the C receiver. Default file: radar_receiver/config/radar_config.bin.
%
%   File format (little-endian): 8 bytes 'RADRXCFG', uint32 version (1),
%   uint32 number of records, then per record: char name[32] (zero padded),
%   uint32 type (0 = real, 1 = complex), uint32 rows, uint32 cols, and
%   rows*cols doubles in column-major order (complex: re, im pairs).

root = fileparts(mfilename('fullpath'));
addpath(root, fullfile(root, 'src'));
if nargin < 1, cfg = radar_config(); end
if nargin < 2, fileName = fullfile(root, 'radar_receiver', 'config', 'radar_config.bin'); end
sys = setup_radar(cfg);

r = {};                                              % {name, value} pairs
r{end+1} = {'nCh', sys.nCh};
r{end+1} = {'nBeams', sys.nBeams};
r{end+1} = {'dbfEnabled', double(sys.dbf.enabled)};
if sys.dbf.enabled
    r{end+1} = {'W', sys.dbf.W};
end
r{end+1} = {'beamU', sys.dbf.beamU};
r{end+1} = {'tilt_deg', sys.ant.tilt_deg};
r{end+1} = {'height_m', sys.ant.height_m};
r{end+1} = {'rate_dps', sys.ant.rate_dps};
r{end+1} = {'azBeamwidth_deg', sys.ant.azBeamwidth_deg};
r{end+1} = {'maxRange_m', cfg.rx.maxRange_m};
r{end+1} = {'nDdc', sys.nDdc};
r{end+1} = {'fsOut_Hz', sys.fsOut};
names = {'short', 'long'};
for i = 1:2
    wv = sys.wave.(names{i});
    r{end+1} = {[names{i} '_ddcFirst'], wv.ddcIdx(1) - 1};      % 0-based
    r{end+1} = {[names{i} '_ddcCount'], numel(wv.ddcIdx)};
    r{end+1} = {[names{i} '_nTaps'], numel(wv.mfTaps)};
    r{end+1} = {[names{i} '_nRange'], wv.nRange};
    r{end+1} = {[names{i} '_taps'], complex(wv.mfTaps)};
    r{end+1} = {[names{i} '_range_m'], wv.range_m};
    r{end+1} = {[names{i} '_delay_s'], wv.delay_s};
    r{end+1} = {[names{i} '_lambda_m'], wv.lambda_m};
end
r{end+1} = {'cpiLength', sys.cpi.length};
r{end+1} = {'nStagger', sys.nStagger};
r{end+1} = {'staggered', double(sys.staggered)};
r{end+1} = {'mtiCoeffs', sys.mtiCoeffs};
r{end+1} = {'nRows', sys.cpi.nRows};
r{end+1} = {'rowLen', sys.cpi.rowLen};
r{end+1} = {'nfft', sys.dop.nfft};
r{end+1} = {'dopWindow', sys.dop.window};
r{end+1} = {'noiseGain', reshape(sys.dop.noiseGain, 1, [])};
r{end+1} = {'dopFreq_Hz', sys.dop.f_Hz};
r{end+1} = {'notchBins', sys.dop.notchBins};
r{end+1} = {'medianToMean', sys.dop.medianToMean};
r{end+1} = {'guardRange', cfg.cfar.guardRange};
r{end+1} = {'trainRange', cfg.cfar.trainRange};
r{end+1} = {'guardDoppler', cfg.cfar.guardDoppler};
r{end+1} = {'trainDoppler', cfg.cfar.trainDoppler};
r{end+1} = {'nTrain', sys.cfar.nTrain};
r{end+1} = {'alpha', sys.cfar.alpha};
r{end+1} = {'sidelobeMargin_dB', cfg.prune.sidelobeMargin_dB};
r{end+1} = {'minHits', max(1, ceil(cfg.extract.minExtent_deg / ...
                       (sys.ant.rate_dps * mean(sys.pri) * cfg.proc.updateEvery)))};
r{end+1} = {'maxMiss', cfg.extract.maxMiss};
r{end+1} = {'updateEvery', cfg.proc.updateEvery};
r{end+1} = {'minExtent_deg', cfg.extract.minExtent_deg};
r{end+1} = {'meanPri_s', mean(sys.pri)};

folder = fileparts(fileName);
if ~isempty(folder) && ~exist(folder, 'dir'), mkdir(folder); end
fid = fopen(fileName, 'w', 'ieee-le');
if fid < 0, error('Cannot write %s', fileName); end
fwrite(fid, 'RADRXCFG', 'char');
fwrite(fid, [1 numel(r)], 'uint32');
for i = 1:numel(r)
    name = r{i}{1};
    value = double(r{i}{2});
    nameBytes = zeros(1, 32);
    nameBytes(1:numel(name)) = double(name);
    fwrite(fid, nameBytes, 'uint8');
    fwrite(fid, [~isreal(value), size(value, 1), size(value, 2)], 'uint32');
    if isreal(value)
        fwrite(fid, value(:), 'double');
    else
        fwrite(fid, [real(value(:)).'; imag(value(:)).'], 'double');
    end
end
fclose(fid);
fprintf('Receiver configuration (%d records, %s PRF) written to %s\n', numel(r), cfg.prf.mode, fileName);
end
