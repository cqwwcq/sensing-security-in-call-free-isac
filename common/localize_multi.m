function mh = localize_multi(cfg, geo, R)
%LOCALIZE_MULTI Multi-hypothesis reconstruction attack (Reviewer 2, Comment 5).
%
%   R is Mtx x Mtx x Ntx, the sample covariances R_{xhat_k} of (3).  Instead
%   of the dominant peak only, the adversary keeps the two strongest distinct
%   peaks of every B_k(theta) (at least one Rayleigh width apart in sin(theta),
%   so a mirror image does not count as a second peak), casts one ray per
%   peak, and takes as its two hypotheses the most-voted cell and the
%   most-voted cell outside a cfg.mh_exclusion-cell neighborhood of it.  The
%   two hypotheses are then resolved by three rules:
%
%   mh.det_votes  more votes wins (ties at random)
%   mh.det_rand   one of the two at random (a mechanism-unaware adversary
%                 that knows one hypothesis may be artificial)
%   mh.det_aware  the hypothesis whose supporting peaks are weaker wins (an
%                 adversary that knows the fake peak is designed to be the
%                 stronger one)
%   mh.det_set    the target cell is one of the two hypotheses, i.e. the
%                 success of an adversary allowed two guesses, which bounds
%                 any rule that picks one of them
%
%   Each det_* is true when the chosen cell is the target's cell.

Ntx = cfg.Ntx; Mtx = cfg.Mtx; cs = cfg.cell_size; A = cfg.area;
theta = linspace(-pi, pi, cfg.n_ang);
At = reshape(beamsteering(theta, Mtx), cfg.n_ang, Mtx);
nx = round((A(2) - A(1)) / cs);
ny = round((A(4) - A(3)) / cs);

votes = zeros(nx, ny);
rays = struct('cells', {}, 'power', {});
for k = 1:Ntx
    B = real(sum((conj(At) * R(:, :, k)) .* At, 2));
    pk = two_peaks(B, theta, Mtx);
    for p = pk
        th = theta(p);
        cand = [th, angle(exp(1j*(pi - th)))];
        valid = false(1, 2);
        cells = cell(1, 2);
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
            rays(end+1) = struct('cells', cells{c}, 'power', B(p)); %#ok<AGROW>
        end
    end
end

% Two hypotheses: the most-voted cell, and the most-voted cell outside its
% neighborhood.
h1 = argmax_random(votes(:));
[i1, j1] = ind2sub([nx, ny], h1);
[I, J] = ndgrid(1:nx, 1:ny);
far = max(abs(I - i1), abs(J - j1)) > cfg.mh_exclusion;
v2 = votes;
v2(~far) = -inf;
h2 = argmax_random(v2(:));
hyp = [h1, h2];

% Mean peak power of the rays supporting each hypothesis.
score = zeros(1, 2);
for h = 1:2
    sup = arrayfun(@(r) any(r.cells == hyp(h)), rays);
    score(h) = mean([rays(sup).power]);
end

tx = min(floor((cfg.PosTarget(1) - A(1)) / cs) + 1, nx);
ty = min(floor((cfg.PosTarget(2) - A(3)) / cs) + 1, ny);
target = sub2ind([nx, ny], tx, ty);

if votes(h1) > votes(h2)
    pick_votes = h1;
else
    pick_votes = hyp(randi(2));
end
if score(1) ~= score(2)
    [~, w] = min(score);
    pick_aware = hyp(w);
else
    pick_aware = hyp(randi(2));
end
mh.det_votes = pick_votes == target;
mh.det_rand = hyp(randi(2)) == target;
mh.det_aware = pick_aware == target;
mh.det_set = any(hyp == target);
end

function pk = two_peaks(B, theta, Mtx)
% Indices of the strongest local maximum and the strongest local maximum at
% least 2/Mtx away from it in sin(theta).
n = numel(B);
prv = B([n, 1:n-1]);
nxt = B([2:n, 1]);
loc = find(B >= prv & B >= nxt);
[~, order] = sort(B(loc), 'descend');
loc = loc(order);
pk = loc(1);
s1 = sin(theta(loc(1)));
for q = loc(2:end).'
    if abs(sin(theta(q)) - s1) >= 2 / Mtx
        pk(2) = q;
        break;
    end
end
end

function i = argmax_random(v)
% Index of the maximum of v, ties broken at random.
best = find(v == max(v));
i = best(randi(numel(best)));
end
