function ref = c_receiver_reference(results)
%C_RECEIVER_REFERENCE  What c_receiver_compare.m needs from the results of run_radar_sim.
%   ref = C_RECEIVER_REFERENCE(results): the MATLAB plots, the track log and
%   the cell sizes (range cell, Doppler bin, FFT length) used to compare.
sys = results.sys;
ref.plots = results.plots;
ref.axes.rangeCell_m = max(diff(sys.wave.long.range_m(1:2)), diff(sys.wave.short.range_m(1:2)));
ref.axes.velocityBin_mps = sys.wave.long.lambda_m / 2 * abs(sys.dop.f_Hz(2) - sys.dop.f_Hz(1));
ref.axes.nfft = sys.dop.nfft;
ref.trackLog = results.tracker.log;
end
