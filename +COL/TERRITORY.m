function dmax_par = TERRITORY(r_par, pp, n_pp)
% "TERRITORY" identifies the aggregates domain boundary as their...
%     ...furthest internal distance with repect to its center of mass.
% ----------------------------------------------------------------------- %
% 
% Inputs:
%     r_par: Center of mass location of independent particles
%     pp: Primary particle information cell array
%     n_pp: Number distribution of primaries
% ----------------------------------------------------------------------- %
% 
% Output:
%     dmax_par: Maximum diameter of independent particles with respect...
%         ...to their center of mass
% ----------------------------------------------------------------------- %

r_par = repelem(r_par, n_pp, 1);
pp = cell2mat(pp);
d_temp = 2 .* sqrt(sum((r_par - pp(:, 3:5)).^2, 2)) + pp(:, 2);...
    % The maximum distance of each primary from its corresponding...
        % ...center of mass
d_temp = mat2cell(d_temp,n_pp);
n_agg = size(n_pp,1); % Total number of aggregates
dmax_par = zeros(n_agg,1);
for i = 1 : n_agg
    dmax_par(i) = max(d_temp{i}); % maximum distance within each...
        % ...aggregate
end

end
