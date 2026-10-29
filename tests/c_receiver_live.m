function results = c_receiver_live(cfg, receiverBin, out)
%C_RECEIVER_LIVE  Run the simulation and stream its DDC output to the C receiver.
%   results = C_RECEIVER_LIVE(cfg, receiverBin, out) exports the receiver
%   configuration of cfg to out/radar_config.bin, starts the C receiver
%   receiverBin in the background on UDP port cfg.udp.port, runs
%   run_radar_sim(cfg) (cfg.udp.enable must be true) and waits until the
%   receiver stops (20 s without data). It returns the results of
%   run_radar_sim. The receiver writes into folder out: ddc_stream.rec (the
%   datagrams, for c_receiver_replay.m), plots_c_live.csv,
%   tracks_c_live.csv and receiver_live.log. Needs a POSIX shell.
if ~exist(out, 'dir'), mkdir(out); end
f = c_receiver_files(out);
for g = {f.rec, f.plotsLive, f.tracksLive, f.liveLog, f.done}
    if exist(g{1}, 'file'), delete(g{1}); end
end
export_receiver_config(cfg, f.cfg);
system(sprintf(['("%s" --config "%s" --port %d --record "%s" --plots-csv "%s" --tracks-csv "%s" ' ...
                '--idle-timeout 20 --quiet > "%s" 2>&1; echo $? > "%s") &'], receiverBin, f.cfg, cfg.udp.port, ...
               f.rec, f.plotsLive, f.tracksLive, f.liveLog, f.done));
pause(2);                                       % the receiver opens its socket
results = run_radar_sim(cfg);
tWait = tic;
while ~exist(f.done, 'file') && toc(tWait) < 120
    pause(1);
end
if ~exist(f.done, 'file')
    error('The C receiver did not stop (see %s)', f.liveLog);
end
fprintf('\n--- C receiver (live) ---\n%s', fileread(f.liveLog));
end
