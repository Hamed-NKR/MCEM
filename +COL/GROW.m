function par = GROW(par)
% "GROW" monitors collision of particles and updates the particle structure
% based on the post-aggregation data.
% ----------------------------------------------------------------------- %

% Input/Output:
    % par: Particle structure
% ----------------------------------------------------------------------- %

n_par = size(par.n,1);  % Total number of particles
dmax = COL.TERRITORY(par.r, par.pp, par.n); % Maximum distance within...
    % ...the aggregates from their center of mass
% Making particle pair indices
ind_pars = (1:n_par)';
ind_pars = [repelem(ind_pars, n_par, 1), repmat(ind_pars, n_par, 1)];

% identifying repeating pairs
rmv1 = (1:n_par)';
rmv1 = repelem((rmv1-1).*n_par, 1:n_par);
rmv2 = repmat((1:n_par)', [1 n_par]);
rmv2 = triu(rmv2);
rmv2 = reshape(rmv2,n_par^2,1);
rmv2(rmv2 == 0) = [];
rmv = rmv1 + rmv2;
ind_pars(rmv,:) = [];

% Generating the "OVR" inputs:
d_pars = [repelem(dmax, n_par, 1), repmat(dmax, n_par, 1)];
    % size input
d_pars(rmv,:) = [];
r_pars = [repelem(par.r, n_par, 1), repmat(par.r, n_par, 1)];
    % location input
r_pars(rmv,:) = [];

% Checking overlapping
ovrs = COL.OVR(r_pars, d_pars);

if ~ isempty(find(ovrs == 1, 1))
    
    % Updating the location of overlapped particles
    ind_chk = ind_pars(ovrs == 1); % Indices of colliding pairs
    n_chk = size(ind_chk,1); % Number of colliding pairs
    
%     % Identifying the repeating elements of the overlapping particle set 
%     [ind_col1, ind_unq1] = unique(ind_chk(:,1), 'stable');
%     ind_dup1 = (setdiff(1:n_chk, ind_unq1))';
%     [ind_col2, ind_unq2] = unique(ind_chk(:,2), 'stable');
%     ind_dup2 = (setdiff(1:n_chk, ind_unq2))';
%     ind_unq = unique(cat(1, ind_col1, ind_col2));
%     [par, ind_col] = COL.UNITE(par, ind_unq); % Merging aggregate pairs...
%         % ...info
    
    ;
end

