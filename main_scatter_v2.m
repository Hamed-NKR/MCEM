clc
clear
close all
warning('off')

%% initialize

% parameters for universal correlation
dpp100 = 17.8;
D_TEM = 0.35;

% parameters for Brasil's correlation
alpha_a = 1.08;
k_a = 1.1;

% limits for plotting Brasil's correlation
npp0_lim = [1e0 2e4];
n_npp0 = 1e4;

% limits for plotting universal correlation
da0_lim = [1e0 2e4];
n_da0 = 1e4;

% scattering parameters
mu = 0; % mean
sigma = 0.2; % sd

%% literature correlations

% dpp vs. npp figure
f1 = figure;
f1.Position = [100, 100, 600, 600];
set(f1, 'color', 'white');

% plot combination of Brasil's and universal correlations
r01 = (npp0_lim(2) / npp0_lim(1)) ^ (1 / (n_npp0 - 1));
npp0 = 1e0 * ones(n_npp0,1) .* r01 .^ (((1 : n_npp0) - 1)');
dpp01 = (((dpp100 ^ (1 / D_TEM)) / 100) * (npp0 / k_a).^(1 / (2 * alpha_a))) .^...
    (D_TEM / (1 - D_TEM));
plt01 = plot(npp0, dpp01, 'Color', [0.4940 0.1840 0.5560],...
    'LineStyle', '-.', 'LineWidth', 2);
hold on

% dpp vs. da figure
f2 = figure;
f2.Position = [100, 100, 600, 600];
set(f2, 'color', 'white');

% plot universal correlation
r02 = (da0_lim(2) / da0_lim(1)) ^ (1 / (n_da0 - 1));
da0 = 1e0 * ones(n_da0,1) .* r02 .^ (((1 : n_da0) - 1)');
dpp02 = dpp100 * (da0 / 100) .^ (D_TEM);
plt02 = plot(da0, dpp02, 'Color', [0.4940 0.1840 0.5560],...
    'LineStyle', '-.', 'LineWidth', 2);
hold on

%% non-scattered library data

% load data from first-stage library
load('D:\Hamed\CND\PhD\My Articles\DLCA1\Results\DAT\database.mat', 'parsdata_sigma')

% plot scaled stage 1 library data in dpp vs. npp domain
figure(f1)
npp_uc = parsdata_sigma{4}(1).npp;
dpp_uc = 1e9 * parsdata_sigma{4}(1).dpp_g(:,1);
plt1_uc = scatter(npp_uc, dpp_uc, 15, hex2rgb('#295F98'), 'o', 'LineWidth', 1);

% plot the same "perfectly" scaled aggregates of above in dpp vs. da domain
figure(f2)
da_uc = 1e9 * parsdata_sigma{4}(1).da;
plt2_uc = scatter(da_uc, dpp_uc, 15, hex2rgb('#295F98'), 'o', 'LineWidth', 1);

%% apply scatter on the perfectly scaled data

% Generate distances that follow a normal distribution
dist_scat = normrnd(mu, sigma, size(npp_uc));

% % parameters of universal correlation in dpp vs. npp domain, in log-log scale
% m = D_TEM / (2 * alpha_a * (1 - D_TEM)); % slope
% b = (((dpp100^(1 / D_TEM)) / 100) * ((1 / k_a)^(1 / (2 * alpha_a)))) ^...
%     (D_TEM / (1 - D_TEM)); % intercept
% 
% % Compute unit vector perpendicular to the universal dpp vs. npp line
% perp_npp = -1 / sqrt(1 + m^2);  % x-component of the perpendicular direction
% perp_dpp = m / sqrt(1 + m^2);   % y-component of the perpendicular direction

% Scatter the dpp vs. npp points vertically from the universal line in...
    % ...dpp vs. npp space
dpp_scat = exp(log(dpp_uc) + dist_scat);

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

figure(f1)

plt1_scat = scatter(pars_scat.n, 1e9 * pars_scat.dpp_g(:,1),...
    15, hex2rgb('#EF9C66'), '^', 'LineWidth', 1);

set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 12,...
    'TickLength', [0.02 0.02], 'XScale', 'log', 'YScale', 'log')
xlim([1,2000])
ylim([5,50])
xlabel('$n_\mathrm{pp}$ [nm]', 'interpreter', 'latex', 'FontSize', 14)
ylabel('$d_\mathrm{pp}$ [nm]', 'interpreter', 'latex', 'FontSize', 14)
legend(cat(2, plt01, plt1_uc, plt1_scat),...
    cat(2, {'Olfert and Rogak (2019)'},...
    {'Original scaling'}, {'Secondary dispersion'}),...
    'interpreter', 'latex', 'FontSize', 12, 'location', 'northwest')

figure(f2)

plt2_scat = scatter(1e9 * pars_scat.da, 1e9 * pars_scat.dpp_g(:,1),...
    15, hex2rgb('#EF9C66'), '^', 'LineWidth', 1);

set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 12,...
    'TickLength', [0.02 0.02], 'XScale', 'log', 'YScale', 'log')
xlim([10,1500])
ylim([5,50])
xlabel('$d_\mathrm{a}$ [nm]', 'interpreter', 'latex', 'FontSize', 14)
ylabel('$d_\mathrm{pp}$ [nm]', 'interpreter', 'latex', 'FontSize', 14)

jj_render = sort(randperm(length(pp_scat), 5));

for j = 1 : length(jj_render)
    figure
    UTILS.PLOTPP(pp_scat{j}(:,3), pp_scat{j}(:,4), pp_scat{j}(:,5),...
        pp_scat{j}(:,2))
end


%% check scaling properties

fit = fitlm(table(log(1e7 * pars_scat.da), log(1e9 * pars_scat.dpp_g(:,1))),...
    'linear', 'Weights', sqrt(pars_scat.n)); % fit a linear regression...
    % ...weighted by sqrt of number of primaries

D_TEM_scat = fit.Coefficients.Estimate(2); % ensemble scaling exponent
dpp100_scat = exp(fit.Coefficients.Estimate(1)); % ensemble scaling prefactor

% 95% confidence intervals for fractal properties
ci = coefCI(fit);
ci_D_TEM_scat = ci(2,:);
ci_dpp100_scat = exp(ci(1,:));

% 95% ci error bars
dci_D_TEM_scat = max(ci_D_TEM_scat) - D_TEM_scat;
dcip_dpp100_scat = max(ci_dpp100_scat) - dpp100_scat;
dcin_dpp100_scat = dpp100_scat - min(ci_dpp100_scat);

% generate the fit data
dpp_scat_fit = dpp100_scat * ((da0/100).^D_TEM_scat);
ci_dpp_scat_fit = [ci_dpp100_scat(1) * (da0/100).^ci_D_TEM_scat(1),...
    ci_dpp100_scat(2) * (da0/100).^ci_D_TEM_scat(2)];

figure(2)
% plot the main fit and CI bounds
plt2_fit = plot(da0, dpp_scat_fit, 'Color', hex2rgb('#C96868'),...
    'LineStyle', '--', 'LineWidth', 1.5);
plt2_err = plot(da0, ci_dpp_scat_fit(:,1), da0, ci_dpp_scat_fit(:,2),...    
    'Color', hex2rgb('#C96868'), 'LineStyle', ':', 'LineWidth', 1);

lgd_fit = strcat('$D_\mathrm{TEM_{scat}}$ =', {' '}, num2str(D_TEM_scat, '%.2f'),...
            {' '}, '$\pm$', {' '}, num2str(dci_D_TEM_scat, '%.2f'), {','},...
            string(newline), '$d_\mathrm{pp,100_{scat}}$ =', {' '}, num2str(dpp100_scat, '%.2f'),...
            {' '}, {'$\pm$'}, {' '}, num2str(max(dcip_dpp100_scat, dcip_dpp100_scat), '%.2f'));

legend(cat(2, plt02, plt2_uc, plt2_scat, plt2_fit),...
    cat(2, {'Olfert and Rogak (2019)'}, {'Original scaling'},...
    {'Secondary dispersion'}, lgd_fit),...
    'interpreter', 'latex', 'FontSize', 12, 'location', 'northwest')

