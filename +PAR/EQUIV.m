function dv = EQUIV(pp)
% "EQUIV" computes the equivalent volumetric diameter of the aggregates.
% ----------------------------------------------------------------------- %
% 
% Inputs:
%     pp: Primary particle information cell array
% ----------------------------------------------------------------------- %
% 
% Output:
%     dv: Volumetric diameter
% ----------------------------------------------------------------------- %

n_agg = size(pp,1); % Total number of aggregates
dv = zeros(n_agg,1); % Initializing the size array

for i = 1 : n_agg
    dv(i) = nthroot(sum(pp{i}(:,2) .^ 3), 3); % Volumetric diameter
end

end
