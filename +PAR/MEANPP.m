function dpp = MEANPP(pp)
% "MEANPP" calculates the mean and standard deviation of primary...
%   ...particle size of aggregates.
% ----------------------------------------------------------------------- %
% 
% Inputs:
%   pp: Primary particle information cell array
% ----------------------------------------------------------------------- %
% 
% Output:
%   dpp: An N*2 array containing:
%       dpp(:,1): Mean primary particle size within each aggregate
%       dpp(:,2): Sample standard deviation
% ----------------------------------------------------------------------- %

n_agg = size(pp,1); % Total number of aggregates
dpp = zeros(n_agg,2); % Initializing the size array

for i = 1 : n_agg
    dpp(i,1) = mean(pp{i}(:,2)); % Average diameter of primaries
    dpp(i,2) = std(pp{i}(:,2)); % Standard deviation
end

end
