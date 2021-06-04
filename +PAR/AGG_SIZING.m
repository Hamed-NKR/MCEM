function d_agg = AGG_SIZING(pp, n_pp)
% "AGG_SIZING" computes various equivalent sizes proposed for fractal...
    % ...aggregates.
% ----------------------------------------------------------------------- %

% Inputs:
    % pp: Properties of primary particle contained within the...
        % ...aggregates
    % n_pp: Number distribution of aggregate primaries
% ----------------------------------------------------------------------- %

% Output:
    % d_agg: Equivalent sizes of the aggregate:
        % Column 1: Mobility diameter
        % Column 2: Aerodynamic diameter
        % Column 3: Equivalent volumetric diameter
% ----------------------------------------------------------------------- %

n_agg = size(n_pp,1); % Total number of aggregates
d_agg = zeros(n_agg,1); % Initializing the aggregates size array

[~, d_agg(:,1)] = COL.EQUIV(pp, n_pp); % Assigning the equivalent...
    % ...volumetric sizes

end

