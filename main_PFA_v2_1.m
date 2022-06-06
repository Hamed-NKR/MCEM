clc
clear
clf('reset')
close all

%% Initializations %%

n_stor = 3; % Number of data storage occurrences
n_try = 3; % Number of DLCA trials

dpp_globstd = 1.4; % Global geometric std of pp size

npp_min = 10; % Aggregate filtering criterion
npp_max = 100; % Iteration limit parameter in terms of number of primaries within the aggregate

j_max = 1e4; % Stage 1 marching index limit

f_dil = 1; % Dilution factor for 2nd stage aggregation

k_max = 1e4; % Stage 2 iteration limit parameter
kk_max = 5; % Stage 2 growth limit parameter


%% 1st stage %%

[params_ud, params_const] = TRANSP.INIT_PARAMS('MCEM_PFAParams'); % Read the input file

pp0 = cell(n_stor, n_try); % Primary particle data storage cell array for initial monodisperse aggregation

pp0_n = cell(n_stor, n_try); % Placeholder for number of primaries withing aggs

i_pp0 = cell(n_stor,1); % Placeholder for index of primary particles

% Iteration steps to be stored
r_str = (npp_max / npp_min)^(1 / (n_stor - 1));
jj = zeros(n_stor,1);
for i = 1 : n_stor
    jj(i) = round(npp_min * (r_str^(i-1)));
end

fprintf('Simulating monodisperse clustering...')
disp(newline)

[pars, fl] = TRANSP.INIT_DOM(params_ud, params_const); % Initialize particle and fluid structs

for i = 1 : n_try
    fprintf('trial %d:', i)
    disp(' ')
    
    % Initilize monodisperse pps for classic DLCA
    [pp_d, pars.n] = PAR.INIT_DIAM(params_ud.Value(5), params_ud.Value(6:7),...
        params_ud.Value(8:10)); % Initialize pp sizes

    pars.pp = mat2cell([(1:size(pp_d))', pp_d, zeros(size(pp_d,1),3),...
        (1:size(pp_d))'], pars.n); % Assign pp indices and sizes

    if params_ud.Value(6) ~= 0
        pars.pp = PAR.INIT_MORPH_RAND(pars.pp); % Randomly initialize pp locations within aggs
    end

    [pars, params_ud] = PAR.INIT_LOC(pars, params_ud); % Randomize initial locations

    pars = PAR.SIZING(pars); % Calculate sizes

    pars = TRANSP.MOBIL(pars, fl, params_const); % Calculate mobility props

    pars.v = PAR.INIT_VEL(pars.pp, pars.n, fl.temp, params_const); % Randomize initial velocities

    j = 2; % initialize iteration index
    jjj = 1; % pp data storage tracking index
    
    fprintf('\b')
    fprintf('DLCA simulation...')
    disp(' ')
    UTILS.TEXTBAR([0, j_max]); % Initialize textbar
    UTILS.TEXTBAR([1, j_max]); % Iteration 1 already done
    
    % 1st stage DLCA iterations
    while (j <= j_max) &&...
            (sum(cat(1, pars.n) >= npp_max) < round(0.95 * length(cat(1, pars.n))))
        [pars, delt] = TRANSP.MARCH(pars, fl, params_const); % Solve for transport
        
        pars = TRANSP.PBC(params_ud.Value(2:4), pars); % Apply periodic BCs
        
        pars = COL.GROW(pars); % Cluster the particles
        
        if sum(cat(1, pars.n) >= jj(jjj)) >= round(0.95 * length(cat(1, pars.n)))
            
            % Store number of primaries
            pp0_n{jjj,i} = cat(1, pars.n);
            
            %         % Remove small aggregates
            %         if jjj == 1
            %             pars.pp(pp0_n{jjj} < npp_min) = [];
            %             pars.n(pp0_n{jjj} < npp_min) = [];
            %             pars.r(pp0_n{jjj} < npp_min, :) = [];
            %             pars.v(pp0_n{jjj} < npp_min, :) = [];
            %             pp0_n{jjj}(pp0_n{jjj} < npp_min) = [];
            %         end
            
            pp0{jjj,i} = pars.pp; % Store pp info
            
            % Update pp indices
            i_pp0{jjj} = cell(length(cat(1, pars.n)), 1);
            for j4 = 1 : length(cat(1, pars.n))
                i_pp0{jjj}{j4} = [pars.pp{j4}(:,1), pars.pp{j4}(:,6)];
                if (i > 1) || (jjj > 1)
                    pp0{jjj,i}{j4}(:,1) = pp0{jjj,i}{j4}(:,1) +...
                        ((i - 1) * n_stor + jjj) * params_ud.Value(5);
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
    if jjj <= n_stor
        pp0{jjj,i} = pars.pp;
        pp0_n{jjj,i} = cat(1, pars.n);
        i_pp0{jjj} = cell(length(cat(1, pars.n)), 1);
        for j4 = 1 : length(cat(1, pars.n))
            i_pp0{jjj}{j4} = [pars.pp{j4}(:,1), pars.pp{j4}(:,6)];
            if (i > 1) || (jjj > 1)
                pp0{jjj,i}{j4}(:,1) = pp0{jjj,i}{j4}(:,1) +...
                    ((i - 1) * n_stor + jjj) * params_ud.Value(5);
            end
        end
    end
    
    % Remove similar aggs
    i_pp0 = cat(1, i_pp0{:}); % Compile indices
    ij = nchoosek(1 : length(cell2mat(pp0_n(:,i))), 2);
    ind_rmv = [];
    for ii = 1 : length(ij)
        if isequal(i_pp0{ij(ii,1)}, i_pp0{ij(ii,2)})
            ind_rmv = [ind_rmv; ij(ii,2)];
        end
    end
    ind_rmv = unique(ind_rmv);
    pp0(ind_rmv) = [];
    pp0_n(ind_rmv) = [];
    
    disp(newline)
end

% Convert the pp data into a 1d cell array
pp0 = pp0(:);
pp0_n = pp0_n(:);

% Remove unused cells
pp0 = pp0(~cellfun('isempty', pp0));
pp0_n = pp0_n(~cellfun('isempty', pp0_n));

pp0 = cat(1, pp0{:}); % Merge pp info from different times
pp0_n = cat(1, pp0_n{:}); % Compile number of primaries data

n_pp0 = sum(pp0_n); % total number of primary particles after stage 1
n_agg0 = size(pp0, 1); % % total number of aggregates before rescaling

dpp0 = PAR.MEANPP(pp0);
dpp0 = dpp0(:,1); % Current mean pp size withing aggs
pp00 = pp0;
dpp_emh = ((17.8^(1/0.35)/100) * (pp0_n / 1.1).^(1 / (2 * 1.08))).^(0.35 / (1 - 0.35)); % Desired mean pp size based on external mixing hypothesis
r_dpp = 1e-9 * dpp_emh ./ dpp0; % Size conversion ratio
% Rescale stored aggregates
for i = 1 : n_agg0
    pp00{i}(:,2:5) = pp00{i}(:,2:5) * r_dpp(i); % Rescale primary particle size
end

pars = structfun(@(x) [], pars, 'UniformOutput', false); % Clear pars struct fields
pars.pp = pp00; % Assign pp info for projected area calculations
pars.n = pp0_n;

disp(' ')
da0 = 2 * sqrt(PAR.PROJECTION(pars, [], 1e4, 20) / pi); % Get projected area diameter for monodisperse populations

dpp00 = PAR.MEANPP(pars.pp);
dpp00 = dpp00(:,1); % Mean primary particle diameter

% Filter particles for a lognormal target area distribution
opts.visual = 'on';
opts.randvar = 'area';
[da1, ind1] = TRANSP.LNSAMPLING(da0, [], [], 10, opts);
da1 = cat(1, da1{:});
ind1 = cat(1, ind1{:});
dpp1 = dpp0(ind1);
pp1 = pp00(ind1);
pp1_n = pp1_n(ind1);

% Randomize order of aggregates
n_agg1 = size(pp1, 1);
i_rnd = randperm(n_agg1);
dpp1 = dpp1(i_rnd);
pp1 = pp1(i_rnd);
pp1_n = pp1_n(i_rnd);

pars.pp = pp1; % Assign pp info for stage 2
pars.n = pp1_n; % Assign number distribution of primaries

%% 2nd stage %%

params_ud.Value(1) = params_ud.Value(1) * f_dil; % Dilute the concentration

[pars, params_ud] = PAR.INIT_LOC(pars, params_ud); % Assign random locations to aggregates

pars = PAR.SIZING(pars); % Get sizes

pars = TRANSP.MOBIL(pars, fl, params_const); % Get mobility props

pars.v = PAR.INIT_VEL(pars.pp, pars.n, fl.temp, params_const); % Assign random velocities to aggregates

opts_grow.indupdate = 'off'; 

disp(newline)
disp('Simulating post-flame mixing...')
UTILS.TEXTBAR([0, k_max]); % Initialize textbar
UTILS.TEXTBAR([1, k_max]); % Iteration 1 already done

k = 2; % initialize iteration index

% Stage 2 DLCA
while (k <= k_max) && (length(cat(1, pars.n)) > round(n_agg0 / kk_max))
    [pars, delt] = TRANSP.MARCH(pars, fl, params_const); % Solve for transport
    
    pars = TRANSP.PBC(params_ud.Value(2:4), pars); % Apply periodic BCs
    
    pars = COL.GROW(pars, opts_grow); % Cluster the particles
    
    % count the number of monodisperse regions within a hybrid...
        %   ...polydisperse aggregates formed by post-flame agglomeration
    pars.n_hyb = COL.HYBRIDITY(pars.pp, pars.n);
    
    pars = PAR.SIZING(pars); % Update sizes
    
    pars = TRANSP.MOBIL(pars, fl, params_const); % Update mobility properties
    
    UTILS.TEXTBAR([k, k_max]); % Update textbar
    
    k = k + 1; % Update iteration ind.
end

da2 = 2 * sqrt(PAR.PROJECTION(pars, [], 1e4, 20) / pi); % Projected area dimater for polydisperse aggs
dpp2 = PAR.MEANPP(pars.pp);
dpp2 = dpp2(:,1);

figure(2)
h2 = gcf;
if ~all(h2.Position == [0, 0, 600, 600])
    h2.Position = [0, 0, 600, 600];
end
set(h2, 'color', 'white');

dpp_uc = linspace(5, 65, 1000);
da_uc = 100 * (dpp_uc / 17.8).^(1 / 0.35);

p21 = plot(da_uc, dpp_uc, 'Color', [0.4660 0.6740 0.1880], 'LineStyle', '-.',...
    'LineWidth', 2.5); % Plot universal correlation
hold on

p22 = scatter(1e9 * da1, 1e9 * dpp1, 25, [0.8500 0.3250 0.0980], 'filled'); % Plot monodisperse aggs

p23 = scatter(1e9 * da2, 1e9 * dpp2, 25, [0 0.4470 0.7410], 'filled'); % Plot hybrid aggs

% axis equal
box on
set(gca, 'FontName', 'SansSerif', 'FontSize', 12, 'TickLength', [0.02 0.02])
xlabel({'\fontsize{14}d_a (nm)', '\fontsize{4} '},'interpreter','tex',...
    'FontName', 'SansSerif', 'FontWeight', 'bold')
ylabel({'\fontsize{4} ', '\fontsize{14}d_p (nm)'},'interpreter','tex',...
    'FontName', 'SansSerif', 'FontWeight', 'bold')
ylim([5, 65])
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
legend([p21, p22, p23], {'Universal correlation', 'Monodisperse', 'Hybrid'},...
    'Location', 'northwest', 'FontName', 'SansSerif', 'FontSize', 12);
title('Primary particle size vs projected area equivalent size',...
    'FontName', 'SansSerif', 'FontWeight', 'bold', 'FontSize', 16)

% figure(3)
% UTILS.RENDER(pars); % display final aggregates

