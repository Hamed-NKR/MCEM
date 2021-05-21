function ind_nn = NNS(par, ind_trg, coef_trg)
% "NNS" identifies the nearest neighbors of a particle based on a...
    % ...user-defined neighboring criterion.

% Inputs:
    % par: The information structure of particles population
    % ind_trg: The index of target particle; i.e., the one neighbors of...
        % ...which need to be identified (an N*1 array)
    % coef_trg: The enlargement coefficient for the size of a spherical...
        % ...barrier used to identify the neighbors
% Output:
    % ind_nn: The index of nearest neighbors

dist_lim = coef_trg .* (par.d(ind_trg)/2) + max(par.d)/2;  % Nearest...
    % ...distance criterion

n_par = size(par.n,1); % Total number of the (independent) particles
n_trg = size(ind_trg, 1); % Number of target particles
ind_chk = (1:n_par)'; % Neighbor checking index list
ind_trg = repmat(ind_trg,n_par,1);

dist_c2c = sqrt(sum((par.r(ind_trg,:) - par.r(ind_chk,:)).^2,2));
    % distance between particles and the target


ind_nn = find(dist_c2c <= dist_lim); % Generating the nearest neighbor list
ind_nn(find(ind_nn == ind_trg(1),1)) = []; % Removing the target...
    % ...particle from the list

end

