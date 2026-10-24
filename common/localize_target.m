function [success, err, est, dom] = localize_target(cfg, geo, R)
%LOCALIZE_TARGET Dominant-peak line fusion of the reconstruction attack.
%
%   R is Mtx x Mtx x Ntx, the sample covariances R_{xhat_k} of (3).  For each
%   AP the dominant angle of B_k(theta) = a^H(theta) R_k a(theta) defines a
%   ray from the AP; the search area is split into square cells, every cell
%   crossed by a ray receives one vote, and the cell with the most votes is
%   the estimate (ties broken at random), as in the earlier multi-static
%   paper.  The ULA cannot tell theta from pi - theta, so only candidates
%   whose ray enters the search area are kept; cfg.mirror_rule = 'random'
%   picks one of them, 'both' casts both rays.
%
%   success  estimate lies in the target's cell
%   err      distance between the estimated cell center and the target
%   dom      per-AP diagnostics: dominant angle, whether the reconstructed
%            fake-direction response exceeds the true-direction response,
%            and where the dominant peak lies (peak_class: 1 fake direction,
%            2 true direction, 3 elsewhere; within half a main lobe,
%            |sin(theta) - sin(phi)| <= 1/Mtx, which also covers the mirror)

Ntx = cfg.Ntx; Mtx = cfg.Mtx; cs = cfg.cell_size; A = cfg.area;
theta = linspace(-pi, pi, cfg.n_ang);
At = reshape(beamsteering(theta, Mtx), cfg.n_ang, Mtx);   % rows a(theta)^T

nx = round((A(2) - A(1)) / cs);
ny = round((A(4) - A(3)) / cs);
votes = zeros(nx, ny);
dom.theta_hat = zeros(1, Ntx);
dom.fake_dominant = false(1, Ntx);
dom.peak_class = zeros(1, Ntx);

for k = 1:Ntx
    Rk = R(:, :, k);
    B = real(sum((conj(At) * Rk) .* At, 2));
    [~, idx] = max(B);
    th = theta(idx);
    if abs(sin(th) - sin(geo.phi(k))) <= 1/Mtx
        dom.peak_class(k) = 2;
    elseif abs(sin(th) - sin(geo.phi_fake(k))) <= 1/Mtx
        dom.peak_class(k) = 1;
    else
        dom.peak_class(k) = 3;
    end
    cand = [th, angle(exp(1j*(pi - th)))];
    cells = cell(1, 2);
    valid = false(1, 2);
    for c = 1:2
        cells{c} = ray_cells(cfg.PosTx(k, :), cand(c), A, cs, nx, ny);
        valid(c) = ~isempty(cells{c});
    end
    if abs(cand(1) - cand(2)) < 1e-9
        valid(2) = false;
    end
    pick = find(valid);
    if numel(pick) > 1 && strcmp(cfg.mirror_rule, 'random')
        pick = pick(randi(numel(pick)));
    end
    for c = pick
        votes(cells{c}) = votes(cells{c}) + 1;
    end
    if ~isempty(pick)
        dom.theta_hat(k) = cand(pick(1));
    end
    bt = real(geo.a_true(:, k)' * Rk * geo.a_true(:, k));
    bf = real(geo.a_fake(:, k)' * Rk * geo.a_fake(:, k));
    dom.fake_dominant(k) = bf > bt;
end

best = find(votes(:) == max(votes(:)));
best = best(randi(numel(best)));
[ix, iy] = ind2sub([nx, ny], best);
est = [A(1) + (ix - 0.5)*cs, A(3) + (iy - 0.5)*cs];

tx = min(floor((cfg.PosTarget(1) - A(1)) / cs) + 1, nx);
ty = min(floor((cfg.PosTarget(2) - A(3)) / cs) + 1, ny);
success = (ix == tx) && (iy == ty);
err = norm(est - cfg.PosTarget);
end

function idx = ray_cells(p0, th, A, cs, nx, ny)
% Linear indices of the cells crossed by the ray from p0 at angle th.
L = hypot(A(2) - A(1), A(4) - A(3));
t = (0:cs/20:L).';
P = p0 + t * [cos(th), sin(th)];
tol = 1e-9;
in = P(:, 1) >= A(1) - tol & P(:, 1) <= A(2) + tol & ...
     P(:, 2) >= A(3) - tol & P(:, 2) <= A(4) + tol;
P = P(in, :);
% A ray from a boundary AP that points outward leaves only the AP itself.
if size(P, 1) < 2
    idx = [];
    return;
end
ix = min(max(floor((P(:, 1) - A(1)) / cs) + 1, 1), nx);
iy = min(max(floor((P(:, 2) - A(3)) / cs) + 1, 1), ny);
idx = unique(sub2ind([nx, ny], ix, iy));
end
