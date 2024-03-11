function [pars, i_list] = UNITE(pars, i_col, col_type)
% "UNITE" unifies the properties of a set of colliding aggregate pairs
% ----------------------------------------------------------------------- %
% 
% Inputs/Outputs:
%     pars: Particle info structure/class
%     i_col: An N*2 (independent) particle index set for colliding pairs
%     i_list: A list of (independent) particle indices before (column 1)...
%         ...and after (column 2) merging
%     col_type: a flag for whether the collision results in...
%         ...aggregation/agglomeration or coalescence
% ----------------------------------------------------------------------- %

% Initial total number of independent particles
if isa(pars, 'AGG')
    n_par0 = length(pars);
else
    n_par0 = size(pars.n, 1);
end

if ~(exist('col_type', 'var')) || isempty(col_type)
    col_type = 'agg';
end

i_list = zeros(n_par0,3); % Initializing the list of indices
i_list(:,1:2) = repmat((1 : n_par0)', 1, 2); % Partcile indices before...
    % ...merging

n_col = size(i_col, 1); % Number of colliding pairs

for i = 1 : n_col
    
    % Initial update of second set of indices
    i_list(find(i_list(:,1) == i_col(i,2), 1), 2) = i_col(i,1);
    
    % Merging primary particle info cells
    if isa(pars, 'AGG')
        pars(i_col(i,1)).pp.r = cat(1, pars(i_col(i,1)).pp.r,...
            pars(i_col(i,2)).pp.r);
        pars(i_col(i,1)).pp.d = cat(1, pars(i_col(i,1)).pp.d,...
            pars(i_col(i,2)).pp.d);
        pars(i_col(i,1)).pp.id = cat(1, pars(i_col(i,1)).pp.id,...
            pars(i_col(i,2)).pp.id);
        
        % Updating the aggregate level properties
        pars(i_col(i,1)).nnl = cat(1, pars(i_col(i,1)).nnl,...
            pars(i_col(i,2)).nnl); % nearest neighbor lists
        pars(i_col(i,1)).nnl((pars(i_col(i,1)).nnl == i_col(i,1)) |...
            (pars(i_col(i,1)).nnl == i_col(i,2))) = []; % Removing the...
                % ...colliding aggregates from neighbor list 
        pars(i_col(i,1)).n = pars(i_col(i,1)).n + pars(i_col(i,2)).n;
            % Number of primaries
        pars(i_col(i,1)).r = AGG.COM(pars(i_col(i,1)).pp); % Center of...
            % ...mass locations
        pars(i_col(i,1)).v = (pars(i_col(i,1)).m .* pars(i_col(i,1)).v +...
            pars(i_col(i,2)).m .* pars(i_col(i,2)).v) ./...
            (pars(i_col(i,1)).m + pars(i_col(i,2)).m); % Net velocties
        pars(i_col(i,1)).m = pars(i_col(i,1)).m + pars(i_col(i,2)).m;
            % Total mass
        
        % Removing the old collided objects
        pars(i_col(i,2)) = [];
        
    else
          pars.pp{i_col(i,1)} = cat(1, pars.pp{i_col(i,1)},...
                pars.pp{i_col(i,2)});
        
%         pars.nnl{i_col(i,1)} = cat(1, pars.nnl{i_col(i,1)},...
%             pars.nnl{i_col(i,2)});
%         pars.nnl{i_col(i,1)}((pars.nnl{i_col(i,1)} == i_col(i,1)) |...
%             (pars.nnl{i_col(i,1)} == i_col(i,2))) = [];
        
        if strcmp(col_type, 'coal') % primary particles merge into a larger sphere (i.e., coalescence)
            
            % equivalent volumetric diameter
            dv = nthroot(sum(pars.pp{i_col(i,1)}(:,2) .^ 3), 3);
            
            % center of mass
            r_com = sum((pars.pp{i_col(i,1)}(:,2) .^ 3) .*...
                pars.pp{i_col(i,1)}(:,3:5), 1) ./...
                sum(pars.pp{i_col(i,1)}(:,2) .^ 3);
            
            pars.pp{i_col(i,1)} = [pars.pp{i_col(i,1)}(1,1), dv, r_com,...
                pars.pp{i_col(i,1)}(1,6)];
            
        end
    end
    
end

i_temp = unique(i_list(:,2), 'stable'); % Extracting remaining old indices
i_temp = [(1 : numel(i_temp))', i_temp]; % Adding new indices
% Asigning new indices to the aggregate populations
for i = 1 : n_par0
    i_list(i,3) = i_temp(find(i_temp(:,2) == i_list(i,2), 1), 1);
end
i_list(:,2) = [];

if ~ isa(pars, 'AGG')
    
    % Updating the particle properties
    if strcmp(col_type, 'coal')
        pars.n(i_col(:,1)) = 1;
        pars.r(i_col(:,1),:) = pars.pp{i_col(:,1)}(1,3:5);
    else
        pars.n(i_col(:,1)) = pars.n(i_col(:,1)) + pars.n(i_col(:,2)); % Number of primaries
        pars.r(i_col(:,1),:) = PAR.COM(pars.pp(i_col(:,1)),...
            pars.n(i_col(:,1))); % Center of mass locations
    end
    
    pars.v(i_col(:,1),:) = (pars.m(i_col(:,1)) .* pars.v(i_col(:,1),:) +...
        pars.m(i_col(:,2)) .* pars.v(i_col(:,2),:)) ./...
        (pars.m(i_col(:,1)) + pars.m(i_col(:,2))); % Net velocties
    pars.m(i_col(:,1)) = pars.m(i_col(:,1)) + pars.m(i_col(:,2)); % Total mass

    % Removing redundant pre-collision data
    pars.pp(i_col(:,2)) = [];
    pars.n(i_col(:,2)) = [];
    pars.r(i_col(:,2),:) = [];
    pars.v(i_col(:,2),:) = [];
    pars.m(i_col(:,2)) = [];
    
    if isfield(pars, 'nnl') && ~isempty(pars.nnl)
        pars.nnl(i_col(:,2)) = [];
    end
    
    % Resetting sizes
    pars.dv = [];
    pars.dg = [];
    pars.dm = [];
    pars.da = [];
    pars.dpp = [];
    pars.dpp_g = [];
    pars.dmax = [];

    % Resetting mass mobility properties
    pars.rho = [];
    pars.delt = [];
    pars.tau = [];
    pars.f = [];
    pars.diff = [];
    pars.lambda = [];
    pars.kn_kin = [];
    pars.kn_diff = [];
    
end

end

