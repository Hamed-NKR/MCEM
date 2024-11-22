clc
clear
clf('reset')
close all

%% Load scaled first-stage LD aggregates %%

% address of data library to be imported
fdir = 'D:\Hamed\CND\PhD\My Articles\DLCA2\mainscatter_sigmapp13';
fname = '18NOV24';
varname = 'pars_out';

% load previously scaled stage 1 aggregate data
load(strcat(fdir, '\', fname, '.mat'), varname)

eval(['pars_LD2' ' = ' varname ';']); % rename loaded structure
 
eval(['clear ', varname]) % remove older name

if ~isfield(pars_LD2, 'pp')
    disp(' ')
    error('Library does not contain aggregates!')
end

%% initialize simulation variables %%

k_max = 1e5; % maximum number of iterations

% assign fractions of aggregates (or times) for second-stage data to be saved
r_n_agg = [1, 0.3, 0.1, 0.03, 0.01];

% resolution of Monte Carlo method for projected area calculation
n_mc_prj = 1e2;
n_ang_prj = 5;

% initial number of aggregates
n0_agg = length(pars_LD2.pp);

% calculate number of primaries in aggregates if missing
if ~isfield(pars_LD2, 'n')

    pars_LD2.n = ones(n0_agg, 1); % allocate space

    for i = 1 : n0_agg
        pars_LD2.n(i) = size(pars_LD2.pp{i}, 1);
    end
end

% calculated projected area diameter if not available
if ~isfield(pars_LD2, 'da')
    pars_LD2.da = 2 * sqrt(PAR.PROJECTION(pars_LD2, [], n_mc_prj,...
        n_ang_prj) / pi);
    disp(' ')
    disp('Calculating projected area...')
end
opts_prj.tbar = 'off';

% calculate characteristic sizes if not input
if ~isfield(pars_LD2, 'dg')
    pars_LD2 = PAR.SIZING(pars_LD2);
end

% read the batch file for fluid and particle input properties
[params_ud, params_const] = TRANSP.INIT_PARAMS('MCEM_PFAParams');

% make the fluid structure
[~, fl] = TRANSP.INIT_DOM(params_ud, params_const);

% change the fluid properties to room condition
opts_fl.amb = 'room';
[fl.mu, fl.lambda] = TRANSP.FLPROPS(fl, params_const, opts_fl);

% calculate initial mobility properties
opts_mobil.mtd = 'continuum'; % choose the method of mobility size calculation
opts_mobil.c_dt = 100;
pars_LD2 = TRANSP.MOBIL(pars_LD2, fl, params_const, opts_mobil);

% adjust the volume fraction 
pp0_ens = cat(1, pars_LD2.pp{:});
params_ud.Value(2) = (pi * sum(pp0_ens(:,2).^3) / (6 * params_ud.Value(1)))^(1/3);
params_ud.Value(3) = params_ud.Value(2);
params_ud.Value(4) = params_ud.Value(2);

% Assign random initial locations and velocities to aggregates
params_ud.Value(1) = 0;
[pars_LD2, params_ud] = PAR.INIT_LOC(pars_LD2, params_ud);
pars_LD2.v = PAR.INIT_VEL(pars_LD2.pp, pars_LD2.n, fl.temp, params_const);

opts_grow.indupdate = 'off'; % flag to update aggregate ids upon each...
    % ...collision (disabled to track hybridty)

% make a placeholder for temporal ensemble data
ensdata.t = zeros(k_max, 1);
ensdata.n_agg = zeros(k_max, 1);
ensdata.tau = zeros(k_max, 2);
ensdata.kn_kin = zeros(k_max, 2);
ensdata.kn_diff = zeros(k_max, 2);
ensdata.dpp = zeros(k_max, 2);
ensdata.sigmapp = zeros(k_max, 2);
% ensdata.da = zeros(k_max, 2);
ensdata.dm = zeros(k_max, 2);

% store ensemble values at the first moment
ensdata.t(1) = 0;
ensdata.n_agg(1) = n0_agg;
ensdata.tau(1,1:2) = [mean(pars_LD2.tau), std(pars_LD2.tau)];
ensdata.kn_kin(1,1:2) = [mean(pars_LD2.kn_kin), std(pars_LD2.kn_kin)];
ensdata.kn_diff(1,1:2) = [mean(pars_LD2.kn_diff), std(pars_LD2.kn_diff)];
ensdata.dpp(1,1:2) = [geomean(pars_LD2.dpp_g(:,1)), UTILS.GEOSTD(pars_LD2.dpp_g(:,1))];
ensdata.sigmapp(1,1:2) = [geomean(pars_LD2.dpp_g(:,2)), UTILS.GEOSTD(pars_LD2.dpp_g(:,2))];
% ensdata.da(1,1:2) = [geomean(pars_LD2.da), UTILS.GEOSTD(pars_LD2.da)];
ensdata.dm(1,1:2) = [geomean(pars_LD2.dm), UTILS.GEOSTD(pars_LD2.dm)];

pars_LD2.n_hyb = ones(n0_agg, 1); % initialize number of sub-aggregates...
    % ...(first, all aggregates are non-hybrid) 

% placeholders for saving second-stage data for individual aggregates
n_dat = length(r_n_agg); % number of moments assigned for data saving
parsdata = struct('dpp', cell(n_dat, 1), 'sigmapp', cell(n_dat, 1),...
    'da', cell(n_dat, 1), 'dg', cell(n_dat, 1), 'npp', cell(n_dat, 1),...
    'n_hyb', cell(n_dat, 1), 'pp', cell(n_dat, 1));

% save data of initial aggregate population
parsdata(1).pp = pars_LD2.pp;
parsdata(1).npp = pars_LD2.n;
parsdata(1).dpp = pars_LD2.dpp_g(:,1);
parsdata(1).sigmapp = pars_LD2.dpp_g(:,2);
% parsdata(1).da = pars_LD2.da;
parsdata(1).dg = pars_LD2.dg;
parsdata(1).n_hyb = pars_LD2.n_hyb;

%% perform Langevin dynamics simulations %%

disp(' ')
disp('Simulating post-flame agglomeration...')
UTILS.TEXTBAR([0, k_max]);
UTILS.TEXTBAR([1, k_max]); % start progress textbar

k = 2; % iteration index

ind_dat = 1; % data saving index

while (k <= k_max) && (ind_dat <= n_dat) && (length(pars_LD2.n) > 1) 
    % check criteria to stop simulations
    
    % solve transport equation
    [pars_LD2, delt] = TRANSP.MARCH(pars_LD2, fl, params_const);
    
    % apply periodic boundary conditions
    pars_LD2 = TRANSP.PBC(params_ud.Value(2:4), pars_LD2);
    
    % join colliding particles
    pars_LD2 = COL.GROW(pars_LD2, opts_grow);
    
    % count number of subaggregates
    pars_LD2.n_hyb = COL.HYBRIDITY(pars_LD2.pp, pars_LD2.n);
    
    % update characteristic sizes
    pars_LD2 = PAR.SIZING(pars_LD2);
    
    % % update projected area sizes
    % pars_LD2.da = 2 * sqrt(PAR.PROJECTION(pars_LD2, [], n_mc_prj,...
    %     n_ang_prj, [], opts_prj) / pi);

    % update mobility properties
    pars_LD2 = TRANSP.MOBIL(pars_LD2, fl, params_const, opts_mobil);
    
    ensdata.t(k) = ensdata.t(k-1) + min(pars_LD2.delt);
    ensdata.n_agg(k) = length(pars_LD2.pp);
    ensdata.tau(k,1:2) = [mean(pars_LD2.tau), std(pars_LD2.tau)];
    ensdata.kn_kin(k,1:2) = [mean(pars_LD2.kn_kin), std(pars_LD2.kn_kin)];
    ensdata.kn_diff(k,1:2) = [mean(pars_LD2.kn_diff), std(pars_LD2.kn_diff)];
    ensdata.dpp(k,1:2) = [geomean(pars_LD2.dpp_g(:,1)), UTILS.GEOSTD(pars_LD2.dpp_g(:,1))];
    ensdata.sigmapp(k,1:2) = [geomean(pars_LD2.dpp_g(:,2)), UTILS.GEOSTD(pars_LD2.dpp_g(:,2))];
    % ensdata.da(k,1:2) = [geomean(pars_LD2.da), UTILS.GEOSTD(pars_LD2.da)];
    ensdata.dm(k,1:2) = [geomean(pars_LD2.dm), UTILS.GEOSTD(pars_LD2.dm)];
    
    if ensdata.n_agg(k) <= (r_n_agg(ind_dat) * ensdata.n_agg(1))
        
        % save data of individual aggregates in selected times
        parsdata(ind_dat).pp = pars_LD2.pp;
        parsdata(ind_dat).npp = pars_LD2.n;
        parsdata(ind_dat).dpp = pars_LD2.dpp_g(:,1);
        parsdata(ind_dat).sigmapp = pars_LD2.dpp_g(:,2);
        % parsdata(ind_dat).da = pars_LD2.da;
        parsdata(ind_dat).dg = pars_LD2.dg;
        parsdata(ind_dat).n_hyb = pars_LD2.n_hyb;
                        
        ind_dat = ind_dat + 1; % update data saving index
        
    end

    if mod(k,250) == 1
        dt = datestr(datetime('now')); % current date and time
        dt = regexprep(dt, ':', '-');
        dt = regexprep(dt, ' ', '_');
        save(strcat('outputs\', dt, '.mat'))
    end

    UTILS.TEXTBAR([k, k_max]); % update progress textbar
    
    k = k + 1; % update iteration index

end

% Remove unused elements from the data storage structures
parsdata(ind_dat:end) = [];
ensdata.t(k:end) = [];
ensdata.n_agg(k:end) = [];
ensdata.tau(k:end,:) = [];
ensdata.kn_kin(k:end,:) = [];
ensdata.kn_diff(k:end,:) = [];
ensdata.dpp(k:end,:) = [];
ensdata.sigmapp(k:end,:) = [];
% ensdata.da(k:end,:) = [];
ensdata.dm(k:end,:) = [];



