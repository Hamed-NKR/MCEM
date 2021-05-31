function [ind_nn, ind_rest] = NNS(par, ind_trg, coef_trg)
% "NNS" identifies the nearest neighbors of a particle based on a...
    % ...user-defined neighboring criterion.

% Inputs:
    % par: The information structure of particles population
    % ind_trg: The index of target particle; i.e., the one neighbors of...
        % ...which need to be identified (an N*1 array)
    % coef_trg: The enlargement coefficient for the size of a spherical...
        % ...barrier used to identify the neighbors (a single value or...
        % ...an n_par*1 array where n_par: number of the particles)
% Outputs:
    % ind_nn: The index of nearest neighbors
    % ind_rest: The index of non-neighbor particles

n_par = size(par.n,1); % Total number of the (independent) particles
n_trg = size(ind_trg,1); % Number of target particles

dmax_par = COL.TERRITORY(par.r, par.pp, par.n);
    % Maximum size of the particles with respect to their center of mass
dist_lim = coef_trg .* dmax_par(ind_trg) ./ 2; % Nearest distance criterion
dist_lim = repelem(dist_lim, n_par, 1);

ind_base = repelem(ind_trg, n_par ,1); % Target checking index list
ind_chk = (1:n_par)'; % Neighbor checking index list
ind_chk = repmat(ind_chk, n_trg, 1);

dist_c2b = sqrt(sum((par.r(ind_base,:) - par.r(ind_chk,:)).^2, 2))...
    - dmax_par(ind_chk); % Distance between the targets center and...
    % ...the other particles boundary

ind_stat = (dist_c2b <= dist_lim) & (ind_chk ~= ind_base);
    % The neighboring status of the particles with respect to the...
        % ...targets (neigbors = 1; non-neighbors & selves = 0)
ind_temp = (dist_c2b > dist_lim) & (ind_chk ~= ind_base);
ind_stat = [ind_chk, ind_stat, ind_temp];
% Adding the particle indices and the non-neighbors status...
    % ...(non-neigbors = 1; neighbors & selves = 0)

% Converting to storable format
ind_stat = mat2cell(ind_stat, n_par .* ones(1,n_trg));
ind_nn = cell(n_trg, 1); % Initializing nearest neighbor list
ind_rest = cell(n_trg, 1); % Initializing non-neighbors list
for i = 1 : n_trg
    ind_nn{i} = ind_stat{i}(ind_stat{i}(:,2) == 1); % Assigning the...
        % ...nearest neighbor indices
    ind_rest{i} = ind_stat{i}(ind_stat{i}(:,3) == 1);
        % Assigning the non-neighbor indices
end

end
