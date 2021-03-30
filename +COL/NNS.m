function ind_nn = NNS(par_r,par_d,ind_trg)
% This function searches for and identifies the nearest neighbors of a
% particle.

% The inputs are the desired particle index and location and size of the...
% particle population.

n_par = size(par_d,1); % Getting number of the particles 
ind_chk = (1:n_par)'; % Neighbor checking index list
ind_chk(ind_trg) = []; % Removing the target particle
ind_trg = repmat(ind_trg,n_par-1,1);

dist = sqrt(sum((par_r(ind_trg,:) - par_r(ind_chk,:)).^2,2));
% distance between particles and the target

par_d(ind_trg) = [];
ind_nn = find(dist <= (20 .* par_d)); % Generating the nearest neighbor list


end

