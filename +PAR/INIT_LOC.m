function pars = INIT_LOC(dom_size, pars)
% "INIT_LOC" randomly distributes the particles throughout the...
%     ...computational domain.
% ----------------------------------------------------------------------- %
% 
% Inputs/Output:
%     dom_size: Computational domain size
%     pars: Particle information structure
% ----------------------------------------------------------------------- %

% Initialization of the location array
n_par = size(pars.n,1); % Total number of (independent) particles
pars.r = PAR.COM(pars.pp, pars.n); % Assigning center of mass as initial...
    % ...location of the aggregates
dmax = PAR.TERRITORY(pars.pp, pars.n); % Maximum distance from the...
    % ...center of each aggregate
pars.r = rand(n_par,3) .* (repmat((dom_size)',n_par,1) -...
     repmat(dmax,1,3))+ (repmat(dmax,1,3) ./ 2);

% Making particle pair indices
ind_pars = (1:n_par)';
ind_pars = [repelem(ind_pars,n_par,1), repmat(ind_pars,n_par,1)];

% identifying repeating pairs
rmv1 = (1:n_par)';
rmv1 = repelem((rmv1-1).*n_par,1:n_par);
rmv2 = repmat((1:n_par)',[1 n_par]);
rmv2 = triu(rmv2);
rmv2 = reshape(rmv2,n_par^2,1);
rmv2(rmv2 == 0) = [];
rmv = rmv1 + rmv2;
ind_pars(rmv,:) = [];

% Generating the "OVR" inputs:
d_pars = [repelem(dmax, n_par, 1), repmat(dmax, n_par, 1)];
    % size input
d_pars(rmv,:) = [];
r_pars = [repelem(pars.r, n_par, 1), repmat(pars.r, n_par, 1)];
    % location input
r_pars(rmv,:) = [];

ovrs = COL.OVR(r_pars, d_pars); % Checking initial overlapping...
% ...between the particles
ind_err = 0; % Initializing error generation index

% Reinitializing overlapped particles
while ~ isempty(find(ovrs == 1, 1))
    
    ind_err = ind_err + 1; % Updating error index
    
    % Updating the location of overlapped particles
    ind_updt = ind_pars(ovrs == 1, 1); % Indices of updated particles
    ind_updt = unique(ind_updt); % removing repeating indices
    pars.r(ind_updt, 1:3) = rand(size(ind_updt,1),3) .*...
        (repmat((dom_size)', size(ind_updt,1),1) -...
        repmat(dmax(ind_updt),1,3)) +...
        (repmat(dmax(ind_updt),1,3) ./ 2);
    r_pars = [repelem(pars.r,n_par,1), repmat(pars.r,n_par,1)];
    r_pars(rmv,:) = [];
    ovrs = COL.OVR(r_pars, d_pars); % Rechecking the overlapping
    
    if ind_err > 10^4
        error('Error assigning random initial locations!\n')
    end
    
end

% Updating the primary particle locations based on their new random...
    % ...center positions
pp_rc = PAR.COM(pars.pp, pars.n); % Initial particle centers
dr0 = pars.r - pp_rc; % Translation vectors
pars.pp = PAR.TRANSLATE(pars.pp, pars.r, pars.n, dr0);

end
