function rmax_agg = TERRITORY(agg_r, agg_pp)
% "TERRITORY" identifies the aggregate domain boundary as its furthest
% internal distance with repect to its center of mass.

% Inputs are the aggregate location and size, and its primary particle
% details.

n_pp = size(agg_pp, 1); % number of primaries within the aggregate
agg_r = repmat(agg_r, n_pp, 1);
rmax_agg = max(sqrt(sum((agg_r - agg_pp(:, 3:6)).^2, 2)) + ...
    agg_pp(:, 2) ./ 2);  % maximum distance of primaries from the...
    % ...aggregate center

end
