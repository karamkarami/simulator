function budgets = print_snr_budget(sys, targets, t_s)
%PRINT_SNR_BUDGET  Print the expected SNR of every target after every processing block.
if nargin < 3, t_s = 0; end
budgets = [];
if isempty(targets), return; end
fprintf('\nExpected SNR per processing block [dB] (CPI centred on the target, t = %.1f s)\n', t_s);
fprintf('%-4s %6s %7s %6s %5s %4s | %7s %7s %6s %6s %6s %8s | %6s %6s\n', 'Tgt', 'RCS', 'R[km]', ...
    'El', 'Pulse', 'Beam', 'ADC', 'DDC', 'DBF', 'MF', 'MTI', 'Doppler', 'Req.', 'Margin');
for i = 1:numel(targets)
    b = snr_budget(sys, targets(i), t_s);
    fprintf('%-4s %6.1f %7.2f %6.2f %5s %4d | %7.1f %7.1f %6.1f %6.1f %6.1f %8.1f | %6.1f %6.1f\n', ...
        b.name, b.rcs_m2, b.range_m/1e3, b.el_deg, b.pulse, b.beam, b.adc, b.ddc, b.dbf, ...
        b.mf, b.mti, b.doppler, b.required, b.margin);
    budgets = [budgets, b]; %#ok<AGROW>
end
fprintf(['ADC: one element, one real sample | DDC: one element, one complex sample | ', ...
         'Doppler: per row, %d row(s) (%.1f independent looks) integrated non-coherently\n'], ...
         sys.cpi.nRows, sys.dop.nLooks);
fprintf('Req.: SNR per row for Pd = 0.9 at Pfa = %.0e (Albersheim, steady target)\n\n', sys.cfar.pfa);
end
