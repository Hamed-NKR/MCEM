clc
clear
clf('reset')
close all

[params_ud, params_const] = TRANSP.INIT_PARAMS('MCEM_DLCAParams'); % Read the input file

n_str = 5; % Number of data storage occurrences
j_max = 1e5; % Marching index limit
jj_max = 200; % Iteration limit parameter

dpp_globstd = 1.4; % global geometric std of pp size

pp0 = cell(n_str,1); % Primary particle data storage cell array for initial monodisperse aggregation

npp_min = 10; % Aggregate filtering criterion
pp0_n = cell(n_str,1); % Placeholder for number of primaries withing aggs

i_pp0 = cell(n_str,1); % Placeholder for index of primary particles

f_dil = 0.01; % Dilution factor for 2nd stage aggregation

% Initilize monodisperse pps for classic DLCA
[pars, fl] = TRANSP.INIT_DOM(params_ud, params_const); % Initialize particle and fluid structs

[pp_d, pars.n] = PAR.INIT_DIAM(params_ud.Value(5), params_ud.Value(6:7),...
    params_ud.Value(8:10)); % Initialize pp sizes

pars.pp = mat2cell([(1:size(pp_d))', pp_d, zeros(size(pp_d,1),3)], pars.n); % Assign pp indices and sizes

if params_ud.Value(6) ~= 0
    pars.pp = PAR.INIT_MORPH_RAND(pars.pp); % Randomly initialize pp locations within aggs
end

[pars, params_ud] = PAR.INIT_LOC(pars, params_ud); % Randomize initial locations

pars = PAR.SIZING(pars); % Calculate sizes

pars = TRANSP.MOBIL(pars, fl, params_const); % Calculate mobility props

pars.v = PAR.INIT_VEL(pars.pp, pars.n, fl.temp, params_const); % Randomize initial velocities

fprintf('Simulating monodisperse clustering...')
disp(' ')
UTILS.TEXTBAR([0, j_max]); % Initialize textbar
UTILS.TEXTBAR([1, j_max]); % Iteration 1 already done

j = 2; % initialize iteration index
jjj = 1; % pp data storage tracking index

% Iterations to be stored
r_str = (jj_max / npp_min)^(1 / (n_str - 1));
jj = zeros(n_str,1);
for i = 1 : n_str
    jj(i) = round(npp_min * (r_str^(i-1))); 
end

% DLCA stage 1 iterations
while (j <= j_max) &&...
        (sum(cat(1, pars.n) >= jj_max) < round(0.95 * length(cat(1, pars.n))))
    [pars, delt] = TRANSP.MARCH(pars, fl, params_const); % Solve for transport
    
    pars = TRANSP.PBC(params_ud.Value(2:4), pars); % Apply periodic BCs
    
    pars = COL.GROW(pars); % Cluster the particles
    
    if sum(cat(1, pars.n) >= jj(jjj)) >= round(0.95 * length(cat(1, pars.n)))
        
        % Store number of primaries
        pp0_n{jjj} = cat(1, pars.n);
        
        % Remove small aggregates
        if jjj == 1
            pars.pp(pp0_n{jjj} < npp_min) = [];
            pars.n(pp0_n{jjj} < npp_min) = [];
            pars.r(pp0_n{jjj} < npp_min, :) = [];
            pars.v(pp0_n{jjj} < npp_min, :) = [];
            pp0_n{jjj}(pp0_n{jjj} < npp_min) = [];
        end
        
        pp0{jjj} = pars.pp; % Store pp info
        
        % Update pp indices
        i_pp0{jjj} = cell(length(cat(1, pars.n)), 1);
        for j4 = 1 : length(cat(1, pars.n))
            i_pp0{jjj}{j4} = pars.pp{j4}(:,1);
            if jjj > 1
                pp0{jjj}{j4}(:,1) = pp0{jjj}{j4}(:,1) + (jjj - 1) * params_ud.Value(5);
            end
        end

        jjj = jjj + 1;
    end
    
    pars = PAR.SIZING(pars); % Update sizes
    
    pars = TRANSP.MOBIL(pars, fl, params_const); % Update mobility properties
    
    UTILS.TEXTBAR([j, j_max]); % Update textbar
    
    j = j + 1; % Update iteration ind.    
end

% Store the pp info and number of pps and update the indices for the last step 
if jjj <= n_str
    pp0{jjj} = pars.pp;
    pp0_n{jjj} = cat(1, pars.n);
    i_pp0{jjj} = cell(length(cat(1, pars.n)), 1);
    for j4 = 1 : length(cat(1, pars.n))
        i_pp0{jjj}{j4} = pars.pp{j4}(:,1);
        if jjj > 1
            pp0{jjj}{j4}(:,1) = pp0{jjj}{j4}(:,1) + (jjj - 1) * params_ud.Value(5);
        end
    end
end

if jjj < n_str
    pp0(jjj + 1 : end) = []; % Remove unused cells
    pp0_n(jjj : end) = [];
    n_str = jjj - 1; % Update the number of stored datasets
end

pp0 = cat(1, pp0{:}); % Merge pp info from different times
pp0_n = cat(1, pp0_n{:}); % Compile number of primaries data
i_pp0 = cat(1, i_pp0{:}); % Compile indices

% Remove similar aggs
ij = nchoosek(1 : length(pp0), 2);
ii_rmv = [];
for i = 1 : length(ij)
    if isequal(i_pp0{ij(i,1)}, i_pp0{ij(i,2)})
        ii_rmv = [ii_rmv; ij(i,2)];
    end
end
ii_rmv = unique(ii_rmv);
pp0(ii_rmv) = [];
pp0_n(ii_rmv) = [];

n_pp0 = sum(pp0_n); % total number of primary particles after stage 1
n_agg0 = size(pp0, 1); % % total number of aggregates before rescaling

pp1 = pp0;
dpp0 = PAR.MEANPP(pp0);
dpp0 = dpp0(:,1); % Current mean pp size withing aggs
dpp_emh = ((17.8^(1/0.35)/100) * (pp0_n / 1.1).^(1 / (2 * 1.08))).^(0.35 / (1 - 0.35)); % Desired mean pp size based on external mixing hypothesis
r_dpp = dpp_emh ./ dpp0; % Size conversion ratio
% Rescale stored aggregates
for i = 1 : n_agg0
    pp1{i}(:,2:5) = pp1{i}(:,2:5) * r_dpp(i); % Rescale primary particle size
end

% Clear pars struct fields
pars = structfun(@(x) [], pars, 'UniformOutput', false);

% Randomize order of aggregates
i_rnd = randperm(n_agg0);
pp1 = pp1(i_rnd);
pp1_n = pp0_n(i_rnd);

pars.pp = pp1; % Assign pp info for stage 2
pars.n = pp1_n; % Assign number distribution of primaries

disp(newline)
da1 = 2 * sqrt(PAR.PROJECTION(pars, [], 1e4, 20) / pi); % Get projected area diameter for monodisperse populations
dpp1 = PAR.MEANPP(pars.pp);
dpp1 = dpp1(:,1); % mean primary particle diameter

params_ud.Value(1) = params_ud.Value(1) * f_dil; % Dilute the concentration
[pars, params_ud] = PAR.INIT_LOC(pars, params_ud); % Assign random locations to aggregates

pars = PAR.SIZING(pars); % Get sizes

pars = TRANSP.MOBIL(pars, fl, params_const); % Get mobility props

pars.v = PAR.INIT_VEL(pars.pp, pars.n, fl.temp, params_const); % Assign random velocities to aggregates

k_max = 1e6; % Iteration limit parameter
kk_max = 5; % Growth limit parameter

disp(newline)
disp('Simulating post-flame mixing...')
UTILS.TEXTBAR([0, k_max]); % Initialize textbar
UTILS.TEXTBAR([1, k_max]); % Iteration 1 already done

k = 2; % initialize iteration index

% Stage 2 DLCA
while (k <= k_max) && (length(cat(1, pars.n)) > round(n_agg0 / kk_max))
    [pars, delt] = TRANSP.MARCH(pars, fl, params_const); % Solve for transport
    
    pars = TRANSP.PBC(params_ud.Value(2:4), pars); % Apply periodic BCs
    
    pars = COL.GROW(pars); % Cluster the particles
    
    pars = PAR.SIZING(pars); % Update sizes
    
    pars = TRANSP.MOBIL(pars, fl, params_const); % Update mobility properties
    
    UTILS.TEXTBAR([k, k_max]); % Update textbar
    
    k = k + 1; % Update iteration ind.
end

da2 = 2 * sqrt(PAR.PROJECTION(pars, [], 1e4, 20) / pi); % Projected area dimater for polydisperse aggs
dpp2 = PAR.MEANPP(pars.pp);
dpp2 = dpp2(:,1);

figure(1)
h2 = gcf;
if ~all(h2.Position == [0, 0, 800, 800])
    h2.Position = [0, 0, 800, 800];
end
set(h2, 'color', 'white');

p21 = scatter(1e9 * dpp1, 1e9 * da1, 25, [0.8500 0.3250 0.0980], 'filled'); % Plot monodisperse aggs
hold on

p22 = scatter(1e9 * dpp2, 1e9 * da2, 25, [0 0.4470 0.7410], 'filled'); % Plot hybrid aggs

dpp_uc = linspace(5, 65, 1000);
da_uc = 100 * (dpp_uc / 17.8).^(1 / 0.35);

p23 = plot(dpp_uc, da_uc, 'Color', [0.4660 0.6740 0.1880], 'LineStyle', '-.',...
    'LineWidth', 2.5); % Plot universal correlation

axis equal
box on
set(gca, 'FontName', 'SansSerif', 'FontSize', 11, 'TickLength', [0.02 0.02])
legend([p21, p22, p23], {'Monodisperse', 'Hybrid', 'Universal correlation'},...
    'FontName', 'SansSerif', 'FontSize', 12);
xlabel({'\fontsize{4} ', '\fontsize{14}d_p (nm)'},'interpreter','tex',...
    'FontName', 'SansSerif', 'FontWeight', 'bold')
ylabel({'\fontsize{14}d_a (nm)', '\fontsize{4} '},'interpreter','tex',...
    'FontName', 'SansSerif', 'FontWeight', 'bold')
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
title('Primary particle size vs projected area equivalent size',...
    'FontName', 'SansSerif', 'FontWeight', 'bold', 'FontSize', 18)

% figure(2)
% UTILS.RENDER(pars); % display final aggregates

