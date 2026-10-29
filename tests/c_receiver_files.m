function f = c_receiver_files(out)
%C_RECEIVER_FILES  Files of a C receiver test in folder out.
f.cfg        = fullfile(out, 'radar_config.bin');     % export_receiver_config.m
f.ref        = fullfile(out, 'matlab_plots.mat');     % MATLAB plots, axes and track log
f.rec        = fullfile(out, 'ddc_stream.rec');       % recorded datagrams
f.plotsLive  = fullfile(out, 'plots_c_live.csv');
f.plotsReplay = fullfile(out, 'plots_c_replay.csv');
f.tracksLive = fullfile(out, 'tracks_c_live.csv');
f.liveLog    = fullfile(out, 'receiver_live.log');
f.replayLog  = fullfile(out, 'receiver_replay.log');
f.done       = fullfile(out, 'receiver_live.done');
end
