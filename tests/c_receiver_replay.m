function logText = c_receiver_replay(receiverBin, cfgFile, recFile, plotsFile, logFile)
%C_RECEIVER_REPLAY  Replay a recorded DDC stream through the C receiver.
%   logText = C_RECEIVER_REPLAY(receiverBin, cfgFile, recFile, plotsFile, logFile)
%   runs the receiver with configuration cfgFile on the datagrams of recFile
%   (c_receiver_live.m), writes its plots to plotsFile and returns its
%   console output (also written to logFile).
if exist(plotsFile, 'file'), delete(plotsFile); end
system(sprintf('"%s" --config "%s" --replay "%s" --plots-csv "%s" --quiet > "%s" 2>&1', ...
               receiverBin, cfgFile, recFile, plotsFile, logFile));
logText = fileread(logFile);
fprintf('\n--- C receiver (replay) ---\n%s', logText);
end
