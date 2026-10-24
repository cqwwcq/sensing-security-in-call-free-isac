function make_detection_results(sweep_file, out_file)
%MAKE_DETECTION_RESULTS Aggregate detection_sweep.mat into detection_results.mat.
%
%   Produces the variables of the published detection_results.mat (p_*,
%   lo_*/hi_* 95% Wilson intervals, eta, eta_zou, cfg) and adds the curves
%   needed for the revision:
%     snr_*_mean   sensing SNR (dB), mean of per-realization dB values as in
%                  the published Fig. 3
%     snrz_prop_mean  proposed design with the AN echo exploited
%     fdom_prop_mean  fraction of APs with a dominant reconstructed fake peak
%     pd_*_mean    mean detection probability, see DETECTION_PD

here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..', 'common'));
if nargin < 1
    sweep_file = fullfile(here, 'detection_sweep.mat');
end
if nargin < 2
    out_file = fullfile(here, 'detection_results.mat');
end
S = load(sweep_file);
cfg = S.cfg;
has_zou = isfield(S, 'det_zou');

O.cfg = cfg;
n = size(S.det_noan, 2);
O.p_noan = mean(S.det_noan);
[O.lo_noan, O.hi_noan] = wilson_ci(O.p_noan, n);
O.p_null = mean(S.det_null, 2).';
[O.lo_null, O.hi_null] = wilson_ci(O.p_null, n);
O.p_prop = mean(S.det_prop, 2).';
[O.lo_prop, O.hi_prop] = wilson_ci(O.p_prop, n);
O.eta = mean(S.eta_prop, 2).';
O.eta_null = mean(S.eta_null, 2).';

O.snr_noan_mean = mean(S.snr_noan);
O.snr_null_mean = mean(S.snr_null, 2).';
O.snr_prop_mean = mean(S.snr_prop, 2).';
O.snrz_prop_mean = mean(S.snrz_prop, 2).';
O.fdom_prop_mean = mean(S.fdom_prop, 2).';
O.pd_noan_mean = mean(S.pd_noan, 3);                 % nOffset x nN
O.pd_null_mean = mean(S.pd_null, 4);                 % nRho x nOffset x nN
O.pd_prop_mean = mean(S.pd_prop, 4);

if has_zou
    nz = size(S.det_zou, 2);
    O.p_zou = mean(S.det_zou, 2).';
    [O.lo_zou, O.hi_zou] = wilson_ci(O.p_zou, nz);
    O.eta_zou = mean(S.eta_zou, 2).';
    O.snr_zou_mean = mean(S.snr_zou, 2).';
    O.pd_zou_mean = mean(S.pd_zou, 4);
end

save(out_file, '-struct', 'O');

% Printout of the values quoted in the paper and the response letter.
[~, j] = min(abs(O.eta - 0.5));
fprintf('No AN: success %.2f %%, sensing SNR %.2f dB\n', 100*O.p_noan, O.snr_noan_mean);
fprintf('eta = %.3f (rho = %.2f):\n', O.eta(j), cfg.rho_list(j));
fprintf('  proposed : success %.2f %%, SNR %.2f dB, SNR with AN echo %.2f dB, fake peak dominant at %.0f %% of APs\n', ...
    100*O.p_prop(j), O.snr_prop_mean(j), O.snrz_prop_mean(j), 100*O.fdom_prop_mean(j));
fprintf('  null-only: success %.2f %%, SNR %.2f dB (its own eta %.3f)\n', ...
    100*O.p_null(j), O.snr_null_mean(j), O.eta_null(j));
if has_zou
    [~, jz] = min(abs(O.eta_zou - 0.5));
    fprintf('  MI-min   : success %.2f %% at eta %.3f, SNR %.2f dB\n', ...
        100*O.p_zou(jz), O.eta_zou(jz), O.snr_zou_mean(jz));
end
i0 = find(cfg.pd_offset_dB == 0, 1);
for b = 1:numel(cfg.pd_N_list)
    fprintf('  Pd (Pfa = %g, N = %d): no AN %.4f, proposed %.4f\n', cfg.Pfa, ...
        cfg.pd_N_list(b), O.pd_noan_mean(i0, b), O.pd_prop_mean(j, i0, b));
end
end
