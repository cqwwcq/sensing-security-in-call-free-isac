function idx = ray_cells(p0, th, A, cs, nx, ny)
%RAY_CELLS Linear indices of the search-area cells crossed by a ray.
%   idx = RAY_CELLS(p0, th, A, cs, nx, ny) samples the ray from p0 at angle
%   th inside the area A = [xmin xmax ymin ymax] with cells of size cs.
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
