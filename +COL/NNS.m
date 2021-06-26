function [ind_nn, varargout] = NNS(par, ind_trg, coef_trg)
% "NNS" identifies the nearest neighbors of a particle based on a...
%     ...user-defined neighboring criterion.
% ----------------------------------------------------------------------- %
% 
% Inputs:
%     par: The information structure of particles population
%     ind_trg: The index of target particle; i.e., the one neighbors of...
%         ...which need to be identified (an N*1 array)
%     coef_trg: The enlargement coefficient for the size of a spherical...
%         ...barrier used to identify the neighbors (a single value or...
%         ...an n_par*1 array where n_par: number of the particles)
% ----------------------------------------------------------------------- %
% 
% Outputs:
%     ind_nn: The index of nearest neighbors
%     varargout{1}: The index of non-neighbor particles
% ----------------------------------------------------------------------- %

if nargout > 2
    error('Error: Invalid number of output arguments!') % Checking for...
        % ...redundant output arguments
end

if isa(par, 'AGG')
    n_par = size(par, 1); % Total number of (independent) particles
    dmax = zeros(n_par);
    for i = 1 : n_par
        dmax(i) = par.TERRITORY(par(i).pp, par(i).r); % Maximum size...
            % ...of the particles with respect to their center of mass
    end
    
else
    n_par = size(par.n, 1);
    dmax = COL.TERRITORY(par.r, par.pp, par.n);
    
end

% Compiling/copying properties locally
r = cat(1, par.r);

n_trg = size(ind_trg,1); % Number of target particles

% Nearest distance criterion
dist_lim = coef_trg .* dmax(ind_trg) ./ 2;
dist_lim = repelem(dist_lim, n_par, 1);

ind_base = repelem(ind_trg, n_par ,1); % Target checking index list
ind_chk = (1:n_par)'; % Neighbor checking index list
ind_chk = repmat(ind_chk, n_trg, 1);

dist_c2b = sqrt(sum((r(ind_base,:) - r(ind_chk,:)).^2, 2))...
    - dmax(ind_chk); % Distance between the targets center and...
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
if nargout > 1
    ind_rest = cell(n_trg, 1); % Initializing non-neighbors list
end
for i = 1 : n_trg
    ind_nn{i} = ind_stat{i}(ind_stat{i}(:,2) == 1); % Assigning the...
        % ...nearest neighbor indices
    if nargout > 1
        ind_rest{i} = ind_stat{i}(ind_stat{i}(:,3) == 1); % Assigning...
            % ...the non-neighbor indices
    end
end

if nargout > 1
    varargout{1} = ind_rest;
end

end
