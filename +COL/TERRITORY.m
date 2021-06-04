function dmax_agg = TERRITORY(r_agg, pp, n_pp)
% "TERRITORY" identifies the aggregate domain boundary as its furthest...
    % ...internal distance with repect to its center of mass.
% ----------------------------------------------------------------------- %

% Inputs:
    % r_agg: Aggregate center of mass location
    % pp: Primary particle information cell array
    % n_pp: Number distribution of primaries
% ----------------------------------------------------------------------- %

% Output:
    % dmax_agg: Maximum aggregate diameter from with respect to its...
        % ...center of mass
% ----------------------------------------------------------------------- %

r_agg = repelem(r_agg, n_pp, 1);
pp = cell2mat(pp);
d_temp = 2 .* sqrt(sum((r_agg - pp(:, 3:5)).^2, 2)) + pp(:, 2);...
    % The maximum distance of each primary from its corresponding...
        % ...center of mass
d_temp = mat2cell(d_temp,n_pp);
n_agg = size(n_pp,1); % Total number of aggregates
dmax_agg = zeros(n_agg,1);
for i = 1 : n_agg
    dmax_agg(i) = max(d_temp{i}); % maximum distance within each...
        % ...aggregate
end

end
