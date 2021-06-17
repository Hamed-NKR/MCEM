function d_par = SIZING(pp, n_pp)
% "AGG_SIZING" computes various equivalent sizes proposed for fractal...
    % ...aggregates.
% ----------------------------------------------------------------------- %

% Inputs:
    % pp: Properties of primary particle contained within the...
        % ...aggregates
    % n_pp: Number distribution of aggregate primaries
% ----------------------------------------------------------------------- %

% Output:
    % d_par: Equivalent sizes of the aggregate:
        % Column 1: Equivalent volumetric diameter
        % Column 2: Gyration diameter
        % Column 3: Mobility diameter
        % Column 4: Aerodynamic diameter
% ----------------------------------------------------------------------- %

n_par = size(n_pp,1); % Total number of aggregates
d_par = zeros(n_par,1); % Initializing the aggregates size array

[~, d_par(:,1)] = COL.EQUIV(pp, n_pp); % Assigning the equivalent...
    % ...volumetric sizes

end

