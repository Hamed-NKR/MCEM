function r_com = COM(pp, n_pp)
% "COM" computes the center of mass of an ensemble of spherical particles.
% ----------------------------------------------------------------------- %
% 
% Inputs:
%     pp: Primary particle information cell array
%     n_pp: Number distribution of primaries
% ----------------------------------------------------------------------- %
% 
% Outputs:
%     r_com: Aggregate center of mass position
% ----------------------------------------------------------------------- %

n_agg = size(n_pp,1); % Total number of aggregates
r_com = zeros(n_agg,3);

for i = 1 : n_agg
    r_com(i,:) = sum((pp{i}(:,2) .^ 3) .* pp{i}(:,3:5), 1) ./...
        sum(pp{i}(:,2) .^ 3); % The center of mass
end

end
