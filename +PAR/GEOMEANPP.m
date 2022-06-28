function dpp_g = GEOMEANPP(pp)
% "GEOMEANPP" calculates the geometric mean and standard deviation of...
%   ...primary particle size of aggregates.
% ----------------------------------------------------------------------- %
% 
% Inputs:
%   pp: Primary particle information cell array
% ----------------------------------------------------------------------- %
% 
% Output:
%   dpp: An N*2 array containing:
%       dpp(:,1): Geometric mean primary particle size within each aggregate
%       dpp(:,2): Sample geometric standard deviation
% ----------------------------------------------------------------------- %

n_agg = size(pp,1); % Total number of aggregates
dpp_g = zeros(n_agg,2); % Initializing the size array

for i = 1 : n_agg
    dpp_g(i,1) = geomean(pp{i}(:,2)); % Geom. average diameter of primaries
    dpp_g(i,2) = UTILS.GEOSTD(pp{i}(:,2)); % Geom. standard deviation
end

end

