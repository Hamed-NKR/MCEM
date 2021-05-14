function ind_nn = NNS(par_r,par_d,ind_trg,coef_trg)
% "NNS" searches for and identifies the nearest neighbors of a particle.

% The inputs are the desired particle index and location and sizes of...
% ...the particle population.

dist_lim = coef_trg * (par_d(ind_trg)/2) + max(par_d)/2;  % Nearest distance limit

n_par = size(par_d,1); % Getting number of the particles 
ind_chk = (1:n_par)'; % Neighbor checking index list
ind_trg = repmat(ind_trg,n_par,1);

dist_c2c = sqrt(sum((par_r(ind_trg,:) - par_r(ind_chk,:)).^2,2));
% distance between particles and the target


ind_nn = find(dist_c2c <= dist_lim); % Generating the nearest neighbor list
ind_nn(find(ind_nn == ind_trg(1),1)) = []; % Removing the target particle...
% ...from the list

end

