function dmax = TERRITORY(pp, n_pp)
% "TERRITORY" identifies the aggregates domain boundary as their...
%     ...furthest internal distance with repect to its center of mass.
% ----------------------------------------------------------------------- %
% 
% Inputs:
%     pp: Primary particle information cell array
%     n_pp: Number distribution of primaries
% ----------------------------------------------------------------------- %
% 
% Output:
%     dmax: Maximum diameter of independent particles with respect...
%         ...to their center of mass
% ----------------------------------------------------------------------- %

n_agg = size(n_pp,1); % Total number of aggregates
dmax = zeros(n_agg,1); % Initializing the size array
r_com = PAR.COM(pp, n_pp); % Center of mass coordinates of aggregates

for i = 1 : n_agg
    dmax(i) = max(2 .* sqrt(sum((repmat(r_com(i,:), n_pp(i), 1) -...
        pp{i}(:,3:5)) .^ 2, 2)) + pp{i}(:,2)); % Maximum distance within...
            % ...each aggregate
end

end
