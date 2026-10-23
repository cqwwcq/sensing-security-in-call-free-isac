function u_new = AN_update_receiver_mmse_null(Ntx, Nue, Mue, Mtx, H_list, w_i, w_t, z, sigma_n, rx_mode) %#ok<INUSL>
%AN_UPDATE_RECEIVER_MMSE_NULL MMSE combiner of eq:mmse under a selectable null set.
%
%   The combiner maximizes the post-combining SINR of user i over
%   {u : Q_i^H u = 0, ||u|| = 1}.  The constraint matrix matches the null
%   constraints that the transmit design actually imposes:
%
%     rx_mode = 0  ->  Q_i = []              (no-AN baseline: no nulls)
%     rx_mode = 1  ->  Q_i = q_iz            (AN-only null, prob:null)
%     rx_mode = 2  ->  Q_i = [q_iz, q_it]    (symmetric null, prob:null and
%                                             prob:null-t)
%
%   Logical input is accepted for backward compatibility: false -> 1,
%   true -> 2.

if nargin < 10 || isempty(rx_mode)
    rx_mode = 2;
end
if islogical(rx_mode)
    rx_mode = 1 + double(rx_mode);
end

u_new = cell(Nue, 1);

% sigma_n may be a scalar, or one effective noise level per user when the
% robust design folds the CSI-error floor L_i^csi into the noise term.
if isscalar(sigma_n)
    sigma_vec = repmat(sigma_n, 1, Nue);
else
    sigma_vec = sigma_n(:).';
end

for i = 1:Nue
    R_i = sigma_vec(i) * eye(Mue);

    for j = 1:Nue
        if j == i
            continue;
        end
        q_ij = zeros(Mue, 1);
        for k = 1:Ntx
            q_ij = q_ij + H_list{i, k} * w_i(:, j, k);
        end
        R_i = R_i + q_ij * q_ij';
    end

    q_it = zeros(Mue, 1);
    q_iz = zeros(Mue, 1);
    q_ii = zeros(Mue, 1);
    for k = 1:Ntx
        q_it = q_it + H_list{i, k} * w_t(:, k);
        q_iz = q_iz + H_list{i, k} * z(:, k);
        q_ii = q_ii + H_list{i, k} * w_i(:, i, k);
    end

    R_i = R_i + q_it * q_it';

    switch rx_mode
        case 0
            Q_i = zeros(Mue, 0);
        case 1
            Q_i = q_iz;
        case 2
            Q_i = [q_iz, q_it];
        otherwise
            error('Unknown rx_mode %d.', rx_mode);
    end

    Rinv_q = R_i \ q_ii;

    if isempty(Q_i)
        u_tilde = Rinv_q;
    else
        Rinv_Q = R_i \ Q_i;
        gram_q = Q_i' * Rinv_Q;
        u_tilde = Rinv_q - Rinv_Q * pinv(gram_q) * (Q_i' * Rinv_q);

        if norm(u_tilde) <= 1e-10 * max(norm(Rinv_q), 1)
            % Numerically stable fallback: solve the MMSE problem directly
            % in the Euclidean null space of Q_i^H.
            N_i = null(Q_i');
            if isempty(N_i)
                error('The null-constrained combiner space is empty for user %d.', i);
            end
            reduced_R = N_i' * R_i * N_i;
            u_tilde = N_i * (reduced_R \ (N_i' * q_ii));
        end
    end

    if norm(u_tilde) <= 0
        u_tilde = ones(Mue, 1) / sqrt(Mue);
    end
    u_new{i} = u_tilde / norm(u_tilde);
end
end
