function [lo, hi] = wilson_ci(p, n, z)
%WILSON_CI Wilson score interval for a binomial proportion p over n trials.
if nargin < 3
    z = 1.96;
end
den = 1 + z.^2 ./ n;
mid = (p + z.^2 ./ (2*n)) ./ den;
half = z .* sqrt(p.*(1 - p)./n + z.^2 ./ (4*n.^2)) ./ den;
lo = max(mid - half, 0);
hi = min(mid + half, 1);
end
