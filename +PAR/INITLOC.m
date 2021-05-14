function pp_r = INITLOC(dom_size,pp_d)
% "INITLOC" randomly distributes the primary particles throughout the
% computational domain.

% Inputs are domain size and primary particle size array (mean+std).

n_pp = size(pp_d,1); % Number of primaries

% Initialization of the location array
pp_r = rand(n_pp,3) .* (repmat((dom_size)',n_pp,1) - repmat(pp_d,1,3)) +...
    (repmat(pp_d,1,3) ./ 2);

% Making particle pair indices
ind_pps = (1:n_pp)';
ind_pps = [repelem(ind_pps,n_pp,1), repmat(ind_pps,n_pp,1)];

% identifying repeating pairs
rmv1 = (1:n_pp)';
rmv1 = repelem((rmv1-1).*n_pp,1:n_pp);
rmv2 = repmat((1:n_pp)',[1 n_pp]);
rmv2 = triu(rmv2);
rmv2 = reshape(rmv2,n_pp^2,1);
rmv2(rmv2 == 0) = [];
rmv = rmv1 + rmv2;
ind_pps(rmv,:) = [];

% Generating the "OVR" inputs:
d_pps = [repelem(pp_d,n_pp,1), repmat(pp_d,n_pp,1)]; % size input
d_pps(rmv,:) = [];
r_pps = [repelem(pp_r,n_pp,1), repmat(pp_r,n_pp,1)]; % location input
r_pps(rmv,:) = [];

ovrs = COL.OVR(r_pps, d_pps); % Checking initial overlapping...
% ...between the primaries
ind_err = 0; % Initializing error generation index

% Reinitializing overlapped particles
while ~ isempty(find(ovrs == 1, 1))
    ind_err = ind_err + 1; % Updating error index
    
    % Updating the location of overlapped pars
    ind_updt = ind_pps(ovrs == 1, 1); % Indices of updated particles
    ind_updt = unique(ind_updt); % removing repeating indices
    pp_r(ind_updt) = rand(size(ind_updt,1),3) .* (repmat((dom_size)',...
        size(ind_updt,1),1) - repmat(pp_d(ind_updt),1,3)) + ...
        (repmat(pp_d(ind_updt),1,3) ./ 2);
    r_pps = [repelem(pp_r,n_pp), repmat(pp_r,n_pp,1)];
    r_pps(rmv,:) = [];
    ovrs = COL.OVR(r_pps, d_pps); % Rechecking the overlapping
    
    if ind_err > 10^2
        error("Error assigning random initial locations!\n")
    end
end

end
