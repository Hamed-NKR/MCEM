function [par, indsout] = UNITE(par, indsin)
% "UNITE" unifies the properties of a set of colliding aggregate pairs
% ----------------------------------------------------------------------- %

% Inputs/Outputs:
    % par: Particle structure
    % indsin: An N*2 particle index set for colliding pairs
    % indsout: The output indices of merged particles
% ----------------------------------------------------------------------- %

n_parin = size(par.n); % Initial total number of independent particles
n_col = size(indsin, 1); % Number of colliding pairs

inds_trc = zeros(n_col,1); % Tracing variable for indexing outputs of...
    % ...collision
for i = 1 : n_col
    inds_trc(i) = par.pp{indsin(i,1)}(1,1);
    par.pp{indsin(i,2)} = COL.CONNECT(par.pp{indsin(i,1)},...
        par.pp{indsin(i,2)});
    par.pp(indsin(i,1)) = [par.pp{indsin(i,1)}; par.pp{indsin(i,2)}];
    
    par.r(indsin(i,1)) = ;
    
end

find

n_pp = par.n(indsin);

for i = 1 : n_col
    pp(indsin) = par.pp{indsin(i,1)};    
end

end

