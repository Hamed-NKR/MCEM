function [pars, i_list] = UNITE(pars, i_col)
% "UNITE" unifies the properties of a set of colliding aggregate pairs
% ----------------------------------------------------------------------- %
% 
% Inputs/Outputs:
%     pars: Particle characteristics structure
%     i_col: An N*2 particle index set for colliding pairs
%     i_list: A list of independent particle indices before (column 1)...
%         ...and after (column 2) merging
% ----------------------------------------------------------------------- %

n_par0 = size(pars.n, 1); % Initial total number of independent particles
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
    pars.pp{i_col(i,1)} = cat(1, pars.pp{i_col(i,1)}, pars.pp{i_col(i,2)});
    
    % Merging nearest neighbor lists
    pars.nnl{i_col(i,1)} = cat(1, pars.nnl{i_col(i,1)}, pars.nnl{i_col(i,2)});
    pars.nnl{i_col(i,1)}((pars.nnl{i_col(i,1)} == i_col(i,1)) |...
        (pars.nnl{i_col(i,1)} == i_col(i,2))) = []; % Removing the...
            % ...colliding aggregates from neighbor list 
end

i_temp = unique(i_list(:,2), 'stable'); % Extracting remaining old indices
i_temp = [(1 : numel(i_temp))', i_temp]; % Adding new indices
% Asigning new indices to the aggregate populations
for i = 1 : n_par0
    i_list(i,3) = i_temp(find(i_temp(:,2) == i_list(i,2), 1), 1);
end
i_list(:,2) = [];

% Updating aggregate level properties
pars.n(i_col(:,1)) = pars.n(i_col(:,1)) + pars.n(i_col(:,2)); % Number...
    % ...of primaries
pars.r(i_col(:,1),:) = PAR.COM(pars.pp(i_col(:,1)), pars.n(i_col(:,1)));
    % Center of mass locations
pars.v(i_col(:,1),:) = (pars.m(i_col(:,1)) .* pars.v(i_col(:,1),:) +...
    pars.m(i_col(:,2)) .* pars.v(i_col(:,2),:)) ./ (pars.m(i_col(:,1)) +...
    pars.m(i_col(:,2))); % Net velocties
pars.m(i_col(:,1)) = pars.m(i_col(:,1)) + pars.m(i_col(:,2)); % Total mass

% Removing redundant pre-collision data
pars.pp(i_col(:,2)) = [];
pars.nnl(i_col(:,2)) = [];
pars.n(i_col(:,2)) = [];
pars.r(i_col(:,2),:) = [];
pars.v(i_col(:,2),:) = [];
pars.m(i_col(:,2)) = [];

% Resetting sizes
pars.dv = [];
pars.dg = [];
pars.dm = [];
pars.da = [];
pars.dpp = [];
pars.dmax = [];

% Resetting mass mobility properties
pars.rho = [];
pars.delt = [];
pars.tau = [];
pars.f = [];
pars.diff = [];
pars.lambda = [];
pars.kn = [];

end

