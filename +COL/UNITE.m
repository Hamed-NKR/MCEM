function [par, i_list] = UNITE(par, i_col)
% "UNITE" unifies the properties of a set of colliding aggregate pairs
% ----------------------------------------------------------------------- %

% Inputs/Outputs:
    % par: Particle characteristics structure
    % i_col: An N*2 particle index set for colliding pairs
    % i_list: A list of independent particle indices before (column 1)...
        % ...and after (column 2) merging
% ----------------------------------------------------------------------- %

n_par0 = size(par.n, 1); % Initial total number of independent particles
i_list = zeros(n_par0,3); % Initializing the list of indices
i_list(:,1:2) = repmat((1 : n_par0)', 1, 2); % Partcile indices before...
    % ...merging
n_col = size(i_col, 1); % Number of colliding pairs
% inds_trc = zeros(n_col, 1); % Introducing a tracing variable for...
%     % ...indexing outputs of collision

for i = 1 : n_col
%     inds_trc(i) = par.pp{i_col(i,1)}(1,1); % Recording the first primary...
%         % ...partcile as the trace
    i_list(find(i_list(:,1) == i_col(i,2), 1), 2) = i_col(i,1);
        % Initial update of second set of indices
    
    % Merging primary particle info cells
    par.pp{i_col(i,1)} = cat(1, par.pp{i_col(i,1)}, par.pp{i_col(i,2)});
    
    % Merging nearest neighbor lists
    par.nnl{i_col(i,1)} = cat(1, par.nnl{i_col(i,1)}, par.nnl{i_col(i,2)});
    par.nnl{i_col(i,1)}((par.nnl{i_col(i,1)} == i_col(i,1)) |...
        (par.nnl{i_col(i,1)} == i_col(i,2))) = []; % Removing the...
            % ...colliding aggregates from neighbor list 
end

i_temp = unique(i_list(:,2), 'stable'); % Extracting remaining old indices
i_temp = [(1 : numel(i_temp))', i_temp]; % Adding new indices
% Asigning new indices to the aggregate populations
for i = 1 : n_par0
    i_list(i,3) = i_temp(find(i_temp(:,1) == i_list(i,2), 1), 2);
end
i_list(:,2) = [];

% Updating aggregate level properties
par.n(i_col(:,1)) = par.n(i_col(:,1)) + par.n(i_col(:,2)); % Number of...
    % ...primaries
[par.r(i_col(:,1),:), par.d(i_col(:,1),1)] = ...
    COL.EQUIV(par.pp(i_col(:,1)), par.n(i_col(:,1))); % Center of mass...
        % ...locations and volumetic sizes
par.v(i_col(:,1),:) = (par.m(i_col(:,1)) .* par.v(i_col(:,1),:) +...
    par.m(i_col(:,2)) .* par.v(i_col(:,2),:)) ./ (par.m(i_col(:,1)) +...
    par.m(i_col(:,2))); % Net velocties
par.m(i_col(:,1)) = par.m(i_col(:,1)) + par.m(i_col(:,2)); % Total mass

% Removing redundant pre-collision data
par.pp(i_col(:,2)) = [];
par.nnl(i_col(:,2)) = [];
par.n(i_col(:,2)) = [];
par.d(i_col(:,2),:) = [];
par.r(i_col(:,2),:) = [];
par.v(i_col(:,2),:) = [];
par.m(i_col(:,2)) = [];

% Resetting mass mobility properties
par.rho = [];
par.delt = [];
par.tau = [];
par.f = [];
par.diff = [];
par.lambda = [];
par.kn = [];

end

