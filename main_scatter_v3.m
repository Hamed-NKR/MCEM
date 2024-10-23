clc
clear
close all
warning('off')

%% initialize, define parameters

% Olfert and Rogak (2019)'s universal correlation
D_TEM = 0.35; % exponent
dpp100 = 17.8; % pefactor
da_lim_uc = [1e0 2e4];  % limits on the projected area diameter (for plotting)
n_da_uc = 1e4; % plot data counts
uc = @(y) dpp100 * (y / 100) .^ D_TEM; % on-demand function for the...
    % ...forward correlation (da and dpp are in [nm])
uc_inv = @(x) 100 * (x / dpp100) .^ (1 / D_TEM); % inverse function to get...
    % ...da from dpp
b0_uc = dpp100 * 100^(-D_TEM);

% Brasil et al. (1999)'s correlation
alpha_a = 1.08; % exponent
k_a = 1.1; % pefactor
npp_lim_bc = [1e0 2e4]; % limits on the number of primaries (for plotting)
n_npp_uc = 1e4; % plot data counts
bc = @(x,y) k_a * (x/y).^(2 * alpha_a); % on-demand function for the...
    % ...forward correlation
bc_inv = @(z) (z / k_a).^(1 / (2 * alpha_a)); % inverse function to get...
    % ...npp from dpp/da

% combination of universal and Brasil correlations
m_ubc = D_TEM / (2 * alpha_a * (1 - D_TEM)); % exponent
b0_ubc = (((dpp100^(1 / D_TEM)) / 100) * ((1 / k_a)^(1 / (2 * alpha_a)))) ^...
    (D_TEM / (1 - D_TEM)); % prefactor
ubc = @(z) b0_ubc * (z .^ m_ubc); % convert, on-demand, number of...
    % ...primaries to primary particle size based on the correlations 
ubc_inv = @(x) (b0_ubc * x) .^ (1 / m_ubc); % inverse combined...
    % ...function to get npp from dpp

% spread variables
mu = 0; % Gaussiam mean of noise distribution around the universal...
    % ...correlation in log-log space
sigma = 0.1; % Gaussian standard deviation
da_noise_lim = [1e1 2e3]; % extents of noise generation 
cn_noise = 10;

% address of aggregate library to be imported for scaling and dispersion
fdir = 'D:\Hamed\CND\PhD\My Articles\DLCA1\Results\DAT\Desktop-simulations\AUG-02-22\sigmapp1\3';
fname = 'wsp_sigmapp_1.3_Aug9';
varname = 'pp0';
vardir = '';

% resolution for projected area calculation
n_mc = 1e2;
n_ang = 5;

%% raw data against the correlations

% load non-scaled first-stage library data
load(strcat(fdir, '\', fname, '.mat'), varname)

% create particle structure  
pars_raw.pp = eval(strcat(varname, vardir)); % store primary particle info
n_agg_raw = length(pars_raw.pp); % number of aggregates
pars_raw.n = zeros(n_agg_raw,1);
for i = 1 : n_agg_raw
    pars_raw.n(i) = size(pars_raw.pp{i}, 1); % count the number of primaries
end
pars_raw = PAR.SIZING(pars_raw); % get characteristic sizes

eval(['clear ', varname])

% initialize scatter data figure
f1 = figure(1);
f1.Position = [50, 50, 800, 450];
set(f1, 'color', 'white');
tt1 = tiledlayout(1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

nexttile(1) % dpp vs. npp figure

% plot combination of Brasil's and universal correlations
r_bc = (npp_lim_bc(2) / npp_lim_bc(1)) ^ (1 / (n_npp_uc - 1));
npp_bc = 1e0 * ones(n_npp_uc,1) .* r_bc .^ (((1 : n_npp_uc) - 1)');
dpp_bc = ubc(npp_bc);
% dpp_bc = (((dpp100 ^ (1 / D_TEM)) / 100) * (npp_bc / k_a).^(1 / (2 * alpha_a))) .^...
%     (D_TEM / (1 - D_TEM));
plt1a_bc = plot(npp_bc, dpp_bc, 'Color', hex2rgb('#597445'),...
    'LineStyle', '-.', 'LineWidth', 2);
hold on

% plot stage-1 library of aggregates in dpp vs npp domain
plt1a_raw = scatter(pars_raw.n, 1e9 * pars_raw.dpp_g(:,1), 8,...
    hex2rgb('#789DBC'), 'o', 'LineWidth', 1);

% appearance configs for dpp vs npp subplot
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 11,...
    'TickLength', [0.02 0.02], 'XScale', 'log', 'YScale', 'log')
xlim([0.8 * min(pars_raw.n), 1.2 * max(pars_raw.n)])
ylim([1e9 * 0.9 * min(pars_raw.dpp_g(:,1)), 1e9 * 1.1 * max(pars_raw.dpp_g(:,1))])
xlabel('$n_\mathrm{pp}$ [-]', 'interpreter', 'latex', 'FontSize', 14)
ylabel('$d_\mathrm{pp}$ [nm]', 'interpreter', 'latex', 'FontSize', 14)

nexttile(2) % dpp vs. da figure

% plot universal correlation
r_uc = (da_lim_uc(2) / da_lim_uc(1)) ^ (1 / (n_da_uc - 1));
da_uc = 1e0 * ones(n_da_uc,1) .* r_uc .^ (((1 : n_da_uc) - 1)');
dpp_uc = uc(da_uc);
plt1b_uc = plot(da_uc, dpp_uc, 'Color', [0.4940 0.1840 0.5560],...
    'LineStyle', '-.', 'LineWidth', 2);
hold on

% plot library in dpp vs. da domain
pars_raw.da = 2 * sqrt(PAR.PROJECTION(pars_raw, [], n_mc, n_ang) / pi);
plt1b_raw = scatter(1e9 * pars_raw.da, 1e9 * pars_raw.dpp_g(:,1), 8,...
    hex2rgb('#789DBC'), 'o', 'LineWidth', 1);

% appearance configs for dpp vs da subplot
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 11,...
    'TickLength', [0.02 0.02], 'XScale', 'log', 'YScale', 'log')
xlim([1e9 * 0.9 * min(pars_raw.da(:,1)), 1e9 * 1.1 * max(pars_raw.da(:,1))])
ylim([1e9 * 0.9 * min(pars_raw.dpp_g(:,1)), 1e9 * 1.1 * max(pars_raw.dpp_g(:,1))])
xlabel('$d_\mathrm{a}$ [nm]', 'interpreter', 'latex', 'FontSize', 14)
ylabel('$d_\mathrm{pp}$ [nm]', 'interpreter', 'latex', 'FontSize', 14)

% generate legend for the entire tile
lgd1 = legend(cat(2, plt1a_bc, plt1b_uc, plt1a_raw),...
    cat(2, {strcat('Brasil et al. (1999) +', string(newline),...
    'Olfert $\&$ Rogak (2019)')}, {'Olfert $\&$ Rogak (2019)'},...
    {strcat('First-stage simulated aggregates', string(newline),...
    '(non-scaled, non-filtered)')}), 'interpreter', 'latex', 'FontSize', 11,...
    'orientation', 'horizontal', 'NumColumns', 3);
lgd1.Layout.Tile = 'south';


%% Generate Gaussian random noise around the universal correlation

n_noise = cn_noise * length(pars_raw.n); % number of noise points for random selection


% % Compute unit vector perpendicular to the universal dpp vs. npp line
% perp_npp = -1 / sqrt(1 + m^2);  % x-component of the perpendicular direction
% perp_dpp = m / sqrt(1 + m^2);   % y-component of the perpendicular direction

% Generate noise values (dispersion distances) that follow a certain distribution
% dist_scat = normrnd(mu, sigma, size(npp_uc)); % Gaussian
dist_scat = mu + sigma * sqrt(12) * (rand(size(npp_uc)) - 0.5); % uniform
% dist_scat = mu + sigma * randn(size(npp_uc)); % Gaussian

% Scatter the dpp vs. npp points vertically from the universal line in...
% ...dpp vs. npp space
dpp_scat = exp(log(dpp_uc) + dist_scat);
% dpp_scat = exp(log(dpp_uc) + dist_scat / cos(atan(m)));

% rescaling factor for noise around universal correlation
rpp_scat = dpp_scat ./ dpp_uc;

pp_scat = parsdata_sigma{4}(1).pp; % load primary particle data

% introduce noise to the perfectly scales aggregates
for i = 1 : length(pp_scat)
    pp_scat{i}(:,2:5) = repmat(rpp_scat(i), npp_uc(i), 4) .* pp_scat{i}(:,2:5);
end

% make a new aggregate population structure for later calculations
pars_scat.n = npp_uc;
pars_scat.pp = pp_scat;

% recalculate projected area based on the the new scaling
pars_scat.da = 2 * sqrt(PAR.PROJECTION(pars_scat, [], 1e2, 5) / pi);

% calculate GM and GSD of primary particle size
pars_scat.dpp_g = PAR.GEOMEANPP(pars_scat.pp);

%% draw scatterd data

nexttile(1)

plt1_scat = scatter(pars_scat.n, 1e9 * pars_scat.dpp_g(:,1),...
    10, hex2rgb('#EF9C66'), '^', 'LineWidth', 1);

set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 11,...
    'TickLength', [0.02 0.02], 'XScale', 'log', 'YScale', 'log')
xlim([1,2000])
ylim([5,50])
xlabel('$n_\mathrm{pp}$ [nm]', 'interpreter', 'latex', 'FontSize', 14)
ylabel('$d_\mathrm{pp}$ [nm]', 'interpreter', 'latex', 'FontSize', 14)

nexttile(2)

plt2_scat = scatter(1e9 * pars_scat.da, 1e9 * pars_scat.dpp_g(:,1),...
    10, hex2rgb('#EF9C66'), '^', 'LineWidth', 1);

set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 11,...
    'TickLength', [0.02 0.02], 'XScale', 'log', 'YScale', 'log')
xlim([10,1500])
ylim([5,50])
xlabel('$d_\mathrm{a}$ [nm]', 'interpreter', 'latex', 'FontSize', 14)
ylabel('$d_\mathrm{pp}$ [nm]', 'interpreter', 'latex', 'FontSize', 14)

%% characterize the imposed dispersion in dpp vs. npp domain

fit1 = fitlm(table(log(pars_scat.n), log(1e9 * pars_scat.dpp_g(:,1))),...
    'linear', 'Weights', sqrt(pars_scat.n)); % fit a linear regression...
% ...weighted by sqrt of number of primaries

D_beta_scat = fit1.Coefficients.Estimate(2); % ensemble scaling exponent
k_beta_scat = exp(fit1.Coefficients.Estimate(1)); % ensemble scaling prefactor

% 95% confidence intervals for fractal properties
ci1 = coefCI(fit1);
ci_D_beta_scat = ci1(2,:);
ci_k_beta_scat = exp(ci1(1,:));

% 95% ci error bars
dci_D_beta_scat = max(ci_D_beta_scat) - D_beta_scat;
dcip_k_beta_scat = max(ci_k_beta_scat) - k_beta_scat;
dcin_k_beta_scat = k_beta_scat - min(ci_k_beta_scat);

% generate the fit data
dpp_scat_fit1 = k_beta_scat * (npp_bc.^D_beta_scat);
ci_dpp_scat_fit1 = [ci_k_beta_scat(1) * (npp_bc.^ci_D_beta_scat(1)),...
    ci_k_beta_scat(2) * (npp_bc.^ci_D_beta_scat(2))];

nexttile(1)
% plot the main fit and CI bounds
plt1_fit = plot(npp_bc, dpp_scat_fit1, 'Color', hex2rgb('#C96868'),...
    'LineStyle', '--', 'LineWidth', 1.5);
plt1_err = plot(npp_bc, ci_dpp_scat_fit1(:,1), npp_bc, ci_dpp_scat_fit1(:,2),...
    'Color', hex2rgb('#C96868'), 'LineStyle', ':', 'LineWidth', 1);

lgd_uc1 = strcat(string(newline), 'Olfert $\&$ Rogak (2019)', string(newline),...
    'and Brasil et al. (1991)', string(newline), '$D_\mathrm{\beta_{uc}}$ =',...
    {' '}, num2str(m, '%.2f'), ', $k_\mathrm{\beta_{uc}}$ =', {' '}, num2str(b, '%.2f'));

lgd_fit1 = strcat(string(newline), string(newline),...
    '$D_\mathrm{\beta_{scat}}$ =', {' '}, num2str(D_beta_scat, '%.2f'),...
    {' '}, '$\pm$', {' '}, num2str(dci_D_beta_scat, '%.2f'), {','},...
    string(newline), '$k_\mathrm{\beta_{scat}}$ =', {' '}, num2str(k_beta_scat, '%.2f'),...
    {' '}, {'$\pm$'}, {' '}, num2str(max(dcip_k_beta_scat, dcin_k_beta_scat), '%.2f'));

legend(cat(2, plt1_uc, plt1_scat, plt1a_bc, plt1_fit),...
    cat(2, {'Original scaling'}, {'Secondary dispersion'},...
    lgd_uc1, lgd_fit1), 'interpreter', 'latex', 'FontSize', 11,...
    'location', 'southoutside', 'orientation', 'horizontal', 'NumColumns', 2)


%% check scaling properties for dpp vs. da

fit2 = fitlm(table(log(1e7 * pars_scat.da), log(1e9 * pars_scat.dpp_g(:,1))),...
    'linear', 'Weights', sqrt(pars_scat.n)); % fit a linear regression...
% ...weighted by sqrt of number of primaries

D_TEM_scat = fit2.Coefficients.Estimate(2); % ensemble scaling exponent
dpp100_scat = exp(fit2.Coefficients.Estimate(1)); % ensemble scaling prefactor

% 95% confidence intervals for fractal properties
ci2 = coefCI(fit2);
ci_D_TEM_scat = ci2(2,:);
ci_dpp100_scat = exp(ci2(1,:));

% 95% ci error bars
dci_D_TEM_scat = max(ci_D_TEM_scat) - D_TEM_scat;
dcip_dpp100_scat = max(ci_dpp100_scat) - dpp100_scat;
dcin_dpp100_scat = dpp100_scat - min(ci_dpp100_scat);

% generate the fit data
dpp_scat_fit2 = dpp100_scat * ((da_uc/100).^D_TEM_scat);
ci_dpp_scat_fit2 = [ci_dpp100_scat(1) * (da_uc/100).^ci_D_TEM_scat(1),...
    ci_dpp100_scat(2) * (da_uc/100).^ci_D_TEM_scat(2)];

nexttile(2)
% plot the main fit and CI bounds
plt2_fit = plot(da_uc, dpp_scat_fit2, 'Color', hex2rgb('#C96868'),...
    'LineStyle', '--', 'LineWidth', 1.5);
plt2_err = plot(da_uc, ci_dpp_scat_fit2(:,1), da_uc, ci_dpp_scat_fit2(:,2),...
    'Color', hex2rgb('#C96868'), 'LineStyle', ':', 'LineWidth', 1);

lgd_uc2 = strcat(string(newline), string(newline), 'Olfert $\&$ Rogak (2019)',...
    string(newline), '$D_\mathrm{TEM}$ = 0.35, $d_\mathrm{pp,100}$ = 17.8');

lgd_fit2 = strcat(string(newline), string(newline), '$D_\mathrm{TEM_{scat}}$ =',...
    {' '}, num2str(D_TEM_scat, '%.2f'), {' '}, '$\pm$', {' '},...
    num2str(dci_D_TEM_scat, '%.2f'), {','}, string(newline),...
    '$d_\mathrm{pp,100_{scat}}$ =', {' '}, num2str(dpp100_scat, '%.2f'),...
    {' '}, {'$\pm$'}, {' '}, num2str(max(dcip_dpp100_scat, dcin_dpp100_scat), '%.2f'));

legend(cat(2, plt1b_uc, plt2_scat, plt1b_uc, plt2_fit),...
    cat(2, {'Original scaling'}, {'Secondary dispersion'},...
    lgd_uc2, lgd_fit2), 'interpreter', 'latex', 'FontSize', 11,...
    'location', 'southoutside', 'orientation', 'horizontal', 'NumColumns', 2)

%% Draw resulting aggregates

f2 = figure(2);
f2.Position = [100, 100, 1200, 800];
set(f2, 'color', 'white');

tt2 = tiledlayout(2, 3, 'Padding', 'none', 'TileSpacing', 'none');

jj_render = sort(randperm(length(pp_scat), 6));

for j = 1 : length(jj_render)

    % figure(j+2)
    % set(gcf, 'Position', [(2*j+3)*50, (2*j-1)*50, 500, 500], 'color', 'white');

    nexttile
    UTILS.PLOTPP(pp_scat{j}(:,3), pp_scat{j}(:,4), pp_scat{j}(:,5),...
        pp_scat{j}(:,2))

end

%% Export plots

exportgraphics(f1, 'outputs\dispersion.png',...
    'BackgroundColor','none', 'ContentType','vector', 'Resolution', 150)
exportgraphics(f2, 'outputs\render.png',...
    'BackgroundColor','none', 'ContentType','vector', 'Resolution', 150)
