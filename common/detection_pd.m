function pd = detection_pd(cfg, V, offset_dB, N_list)
%DETECTION_PD Detection probability of the GLRT for given transmit beams.
%
%   pd = DETECTION_PD(cfg, V, offset_dB, N_list) returns a numel(offset_dB)
%   x numel(N_list) matrix.  V (Nrx x Ntx x Ns) holds the echo vectors from
%   SENSING_METRICS.  The bistatic coefficients alpha_{r,k} ~ CN(0, 1) are
%   drawn independently per transmit-receive pair (Swerling I, random echo
%   phases), so the per-receive-AP SNR is gamma_r = ||sum_k alpha_{r,k}
%   V(r,k,:)||^2 / sigma_rx^2, whose mean is the per-receive-AP SNR behind (2).
%
%   Each receive AP matched-filters the Ns = Nue+1 known streams, and the
%   GLRT sums the output energies over streams and receive APs.  The
%   normalized statistic is chi-square with 2 Nrx Ns degrees of freedom
%   under H0 and noncentral chi-square with noncentrality 2 N sum_r gamma_r
%   under H1.  offset_dB scales all echoes, e.g. to emulate a smaller RCS.

[Nrx, Ntx, Ns] = size(V);
sigma_rx = 10^(cfg.sigma_rx_dBm/10) / 1000;
df = 2 * Nrx * Ns;
tau = chi2inv(1 - cfg.Pfa, df);

D = cfg.pd_alpha_draws;
gsum = zeros(D, 1);
for r = 1:Nrx
    alpha = (randn(D, Ntx) + 1j*randn(D, Ntx)) / sqrt(2);
    cr = alpha * reshape(V(r, :, :), Ntx, Ns);          % D x Ns
    gsum = gsum + sum(abs(cr).^2, 2) / sigma_rx;
end

pd = zeros(numel(offset_dB), numel(N_list));
for a = 1:numel(offset_dB)
    for b = 1:numel(N_list)
        lam = 2 * N_list(b) * 10^(offset_dB(a)/10) * gsum;
        pd(a, b) = mean(ncx2cdf(tau, df, lam, 'upper'));
    end
end
end
