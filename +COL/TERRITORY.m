function dmax_agg = TERRITORY(agg_r, agg_ppinfo)
% "TERRITORY" identifies the aggregate domain boundary as its furthest
% internal distance with repect to its center of mass.

% Inputs are the aggregate location and size, and its primary particle
% details.

n_pp = size(agg_ppinfo, 1); % number of primaries within the aggregate
agg_r = repmat(agg_r, n_pp, 1);
dist = sqrt(sum((agg_r - agg_ppinfo(:, 3:6)).^2, 2)) + ...
    agg_ppinfo(:, 2) ./ 2;  % maximum distance of primaries from the...
    % ...aggregate center
dmax_agg = 2 * max(dist);

end

