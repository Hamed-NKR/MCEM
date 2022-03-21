clc
clear
clf('reset')
close all

[params_ud, params_const] = TRANSP.INIT_PARAMS('MCEM_PFAParams'); % Read the input file

mixlev = 10; % Number of monodisperse DLCA runs to produce PFA polydispersity
j_max = 2e4 * ones(mixlev,1); % Marching index limit
pp_pfa0 = cell(mixlev,1); % Primary particle storage variable
dpp_std = 1e-8; % pp size std between the iterations
dpp_mean = params_ud.Value(7); % pp size mean btw. ~
del_dpp = randn(mixlev,1); % normally random selected dpp size params
npp0 = params_ud.Value(4); % number of primaries in each population

for i = 1 : mixlev    
    time = zeros(j_max(i),1); % Aggregation time array
    
    params_ud.Value(7) = dpp_std .* del_dpp(i) + dpp_mean; % Assign normal dist. to the pp monodisperse size

    % Initilize monodisperse pps for classic DLCA
    [pars, fl] = TRANSP.INIT_DOM(params_ud, params_const); % Initialize particle and fluid structs
    
    [pp_d, pars.n] = PAR.INIT_DIAM(params_ud.Value(4), params_ud.Value(5:6),...
        params_ud.Value(7:9)); % Initialize pp sizes
    
    pars.pp = mat2cell([(1:size(pp_d))', pp_d, zeros(size(pp_d,1),3)], pars.n); % Assign pp indices and sizes
    
    if params_ud.Value(6) ~= 0
        pars.pp = PAR.INIT_MORPH_RAND(pars.pp); % Randomly initialize pp locations within aggs
    end
    
    pars = PAR.INIT_LOC(params_ud.Value(1:3), pars); % Randomize initial locations
    
    pars = PAR.SIZING(pars); % Calculate sizes
    
    pars = TRANSP.MOBIL(pars, fl, params_const); % Calculate mobility props
    
    pars.v = PAR.INIT_VEL(pars.pp, pars.n, fl.temp, params_const); % Randomize initial velocities
    
    fprintf('Simulating population %d:', i)
    UTILS.TEXTBAR([0, j_max(i)]); % Initialize textbar
    UTILS.TEXTBAR([1, j_max(i)]); % Iteration 1 already done
    
    j = 2; % initialize iteration index
    
    % DLCA stage 1 iterations
    while (j <= j_max(i)) && (all(cat(1, pars.n) < (params_ud.Value(4) / 3)))
        [pars, delt] = TRANSP.MARCH(pars, fl, params_const); % Solve for transport
        
        pars = TRANSP.PBC(params_ud.Value(1:3), pars); % Apply periodic BCs
        
        pars = COL.GROW(pars); % Cluster the particles
        
        pars = PAR.SIZING(pars); % Update sizes
        
        pars = TRANSP.MOBIL(pars, fl, params_const); % Update mobility properties
        
        time(j) = time(j-1) + delt; % Update time
        
        UTILS.TEXTBAR([j, j_max(i)]); % Update textbar
        
        j = j + 1; % Update iteration ind.
    end
    
    if i >= 2
        for ii = 1 : length(pars.n)
            pars.pp{ii}(:,1) = pars.pp{ii}(:,1) + (i-1) * npp0; % update pp indices
        end
    end
    
    pp_pfa0{i} = pars.pp; % store pp info
end

% clear pars struct fields
pars = structfun(@(x) [], pars, 'UniformOutput', false);

% Initialize PFA parameters

l_dom = 2 * (mixlev)^(1/3) * params_ud.Value(1:3); % Update domain size

pp_pfa = []; % Initialize stage 2 agg pp info
for ii = 1 : mixlev
    pp_pfa = [pp_pfa; pp_pfa0{ii}]; % Concatinate stage 1 aggregate pp info
end

n_agg0 = length(pp_pfa); % Initial numer of aggregates
pp_pfa = pp_pfa(randperm(n_agg0)); % Randomize order of aggregates

% Update number of primaries field
pars.n = zeros(n_agg0,1);
for ii = 1 : n_agg0
    pars.n(ii) = size(pp_pfa{ii},1);
end

pars.pp = pp_pfa; % Assign pp info for stage 2

pars = PAR.INIT_LOC(l_dom, pars); % Assign random locations to aggregates

pars = PAR.SIZING(pars); % Get sizes
    
pars = TRANSP.MOBIL(pars, fl, params_const); % Get mobility props

pars.v = PAR.INIT_VEL(pars.pp, pars.n, fl.temp, params_const); % % Assign random velocities to aggregates

k_max = 2e4;
time = zeros(k_max,1); % Initialize time array

disp('Simulating Post-flame Agglomeration:')
UTILS.TEXTBAR([0, k_max]); % Initialize textbar
UTILS.TEXTBAR([1, k_max]); % Iteration 1 already done

k = 2; % initialize iteration index

% Stage 2 DLCA
while (k <= k_max) && (length(cat(1, pars.n)) > (n_agg0 / 5))
    [pars, delt] = TRANSP.MARCH(pars, fl, params_const); % Solve for transport
    
    pars = TRANSP.PBC(l_dom, pars); % Apply periodic BCs
    
    pars = COL.GROW(pars); % Cluster the particles
    
    pars = PAR.SIZING(pars); % Update sizes
    
    pars = TRANSP.MOBIL(pars, fl, params_const); % Update mobility properties
    
    time(k) = time(k-1) + delt; % Update time
    
    UTILS.TEXTBAR([k, k_max]); % Update textbar
    
    k = k + 1; % Update iteration ind.
end

figure(1)
UTILS.RENDER(pars); % display final aggregates


