function cfg = apply_overrides(cfg, overrides)
%APPLY_OVERRIDES Replace fields of cfg by those of the struct overrides.
%   Lets a driver be run with different settings without editing it, e.g.
%   run_detection_sweep('', struct('num_runs', 1, 'use_parfor', false)).
if nargin < 2 || isempty(overrides)
    return;
end
for f = fieldnames(overrides).'
    if ~isfield(cfg, f{1})
        error('apply_overrides: unknown setting "%s".', f{1});
    end
    cfg.(f{1}) = overrides.(f{1});
end
end
