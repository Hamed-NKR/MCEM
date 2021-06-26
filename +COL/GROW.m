function par = GROW(par)
% "GROW" monitors collision of particles and updates the particle...
%     ...structure based on the post-aggregation data.
% ----------------------------------------------------------------------- %
% 
% Input/Output:
%     par: Particle structure
% ----------------------------------------------------------------------- %

n_par = size(par.n,1);  % Total number of particles
dmax = COL.TERRITORY(par.r, par.pp, par.n); % Maximum distance within...
    % ...the aggregates from their center of mass
% Making particle pair indices
ind_pars = (1 : n_par)';
ind_pars = [repelem(ind_pars, n_par, 1), repmat(ind_pars, n_par, 1)];

% identifying repeating pairs
rmv1 = (1 : n_par)';
rmv1 = repelem((rmv1 - 1) .* n_par, 1 : n_par);
rmv2 = repmat((1 : n_par)', [1 n_par]);
rmv2 = triu(rmv2);
rmv2 = reshape(rmv2, n_par^2, 1);
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

% Updating the location of overlapped particles
if ~ isempty(find(ovrs == 1, 1))
    
    ind_chk = ind_pars(ovrs == 1, :); % Indices of colliding pairs
    n_chk = size(ind_chk,1); % Number of colliding pairs
    
    for i = 1 : n_chk
        if ind_chk(i,1) ~= ind_chk(i,2)

            % Sticking the two colliding particles
            [par.pp{ind_chk(i,1)}, par.pp{ind_chk(i,2)}] =...
                COL.CONNECT_SIMPLE(par.pp{ind_chk(i,1)},...
                par.pp{ind_chk(i,2)}); 
%             [par.pp{ind_chk(i,1)}, par.pp{ind_chk(i,2)}, colstat] =...
%                 COL.CONNECT_REL(par.pp{ind_chk(i,1)},...
%                 par.v(ind_chk(i,1),:), par.pp{ind_chk(i,2)},...
%                 par.v(ind_chk(i,2),:)); 
%             par.pp{ind_chk(i,2)} = COL.CONNECT_RAND(par.pp{ind_chk(i,1)}...
%                 , par.pp{ind_chk(i,2)}); 
            
            % Defining collision status variable(1 --> collided. 0 -->...
                % ...uncollided)
            if ~exist('colstat', 'var'); colstat = 1; end
            
            if colstat
                
                % Merging the particles info
                [par, ind_new] = COL.UNITE(par, [ind_chk(i,1),...
                    ind_chk(i,2)]);

                % Updating the general overlap checking indices after each...
                    % ...unification
                for j = 1 : n_chk
                    ind_chk(j,1) = ind_new(find(ind_new(:,1) ==...
                        ind_chk(j,1), 1), 2);
                    ind_chk(j,2) = ind_new(find(ind_new(:,1) ==...
                        ind_chk(j,2), 1), 2);
                end
            end
        end
    end
    
end

