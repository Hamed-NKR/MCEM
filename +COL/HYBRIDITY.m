function n_hyb = HYBRIDITY(pp, n_pp)
% "HYBRIDITY" counts the number of monodisperse regions within a...
%   ...hybrid plydisperse aggregates formed by post-flame agglomeration.
% ----------------------------------------------------------------------- %
% 
% Inputs:
%     pp: Primary particle information cell array
%     n_pp: Number distribution of primaries
% ----------------------------------------------------------------------- %
% 
% Output:
%     n_hyb: number of distincts regions in hybrid aggregates
% ----------------------------------------------------------------------- %

n_hyb = zeros(n_agg,1); % Initializing the number of regions array

for i = 1 : length(n_pp)
    n_hyb(i) = length(unique(pp{i}(:,6)));
end

end

