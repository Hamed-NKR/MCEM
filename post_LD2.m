clc
clear
clf('reset')
close all

%% Load second-stage LD aggregates %%

% address of second-stage data to be imported
fdir = 'D:\Hamed\CND\PhD\My Articles\DLCA2\mainscatter_sigmapp10\FLAT';
fname = 'LD2-26-Nov-2024_23-56-35_Final';

% variables of interest in the data
varnames = {'parsdata', 'ensdata', 'r_n_agg', 'fl'}; 

% load simulated data
for i = 1 : numel(varnames)
    load(strcat(fdir, '\', fname, '.mat'), varnames{i})
end

n_dat = length(r_n_agg); % number of data times to be plotted

% initialize temporal dpp vs da figure
f1 = figure(1);
f1.Position = [50, 50, 500, 600];
set(f1, 'color', 'white')

% placholders for plots & legends
plt1 = cell(n_dat+1, 1);
legtxt1 = cell(n_dat+1, 1);
legtxt1{end} = 'Olfert $\&$ Rogak (2019)';

% initialize marker colors, shapes and sizes
mc1 = colormap(hot);
cind1 = round(1 + (length(mc1) - 1) .* (0.05 : 0.7 / (n_dat - 1) : 0.75)');
mc1 = mc1(cind1,:);
ms1 = 10 * ones(n_dat, 1);
ms1(4) = 15;
ms1(5) = 12;
mt1 = {'o', '^', 'v', 's', 'd', 'p'};
mc1(end,:) = [236,230,61] / 255;

% generate and plot universal correlation for dpp vs. da
D_TEM = 0.35; % exponent
dpp_100 = 17.8; % pefactor
da_lim_uc = [1e0 2e4];  % limits on the projected area diameter
n_da_uc = 1e4; % number of data
uc1 = @(y) dpp_100 * (y / 100) .^ D_TEM; % on-demand function for the...
    % ...forward correlation in the geometrical domain (dpp as a...
    % ...function of da in [nm])
r_uc1 = (da_lim_uc(2) / da_lim_uc(1)) ^ (1 / (n_da_uc - 1));
da_uc = da_lim_uc(1) * ones(n_da_uc,1) .* r_uc1 .^ (((1 : n_da_uc) - 1)');
dpp_uc = uc1(da_uc);
plt1{end} = plot(da_uc, dpp_uc, 'Color', [0.4940 0.1840 0.5560],...
    'LineStyle', '-.', 'LineWidth', 3);
hold on

for i = 1 : n_dat
    
    % plot temporal data for dpp vs. da
    plt1{i} = scatter(1e9 * parsdata(i).da, 1e9 * parsdata(i).dpp,...
        ms1(i), mc1(i,:), mt1{i}, 'LineWidth', 1); 
    
    % make legends and adjust their format
    if i == 1
        legtxt1(i) = strcat('$n_\mathrm{agg}/n_\mathrm{agg_0}$ =',...
            {' '}, num2str(r_n_agg(i), '%.0f'));
    elseif ismember(i, [2,3])
        legtxt1(i) = strcat('$n_\mathrm{agg}/n_\mathrm{agg_0}$ =',...
            {' '}, num2str(r_n_agg(i), '%.1f'));
    else
        legtxt1(i) = strcat('$n_\mathrm{agg}/n_\mathrm{agg_0}$ =',...
            {' '}, num2str(r_n_agg(i), '%.2f'));
    end
end

dx1 = [1e9 * 0.9 * min(cat(1, parsdata.da)),...
    1e9 * 1.1 * max(cat(1, parsdata.da))];
dy1 = [1e9 * 0.9 * min(cat(1, parsdata.dpp)),...
    1e9 * 1.1 * max(cat(1, parsdata.dpp))];

% plot appearance configs
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 11,...
    'TickLength', [0.02 0.02], 'XScale', 'log', 'YScale', 'log')
xlim(dx1)
ylim(dy1)
xlabel('$d_\mathrm{a}$ [nm]', 'interpreter', 'latex', 'FontSize', 14)
ylabel('$d_\mathrm{pp}$ [nm]', 'interpreter', 'latex', 'FontSize', 14)
box on

% generate legends
lgd1 = legend(cat(1, plt1{:}), legtxt1, 'interpreter', 'latex',...
    'FontSize', 11, 'orientation', 'horizontal', 'NumColumns', 2,...
    'Location', 'southoutside');


% initialize temporal rho vs dm figure
f2 = figure(2);
f2.Position = [100, 100, 500, 600];
set(f2, 'color', 'white')

% placholders for plots & legends
plt2 = cell(n_dat+1, 1);
legtxt2 = cell(n_dat+1, 1);
legtxt2{end} = 'Olfert $\&$ Rogak (2019)';

% initialize marker colors, shapes and sizes
mc2 = colormap(hot);
cind2 = round(1 + (length(mc2) - 1) .* (0.05 : 0.7 / (n_dat - 1) : 0.75)');
mc2 = mc2(cind2,:);
ms2 = 10 * ones(n_dat, 1);
ms2(4) = 15;
ms2(5) = 12;
mt2 = {'o', '^', 'v', 's', 'd', 'p'};
mc2(end,:) = [236,230,61] / 255;

% generate and plot universal correlation for rho_eff vs. dm
D_m = 2.48; % exponent
rho_eff_100 = 510; % pefactor
dm_lim_uc = [1e0 2e4];  % limits on the mobility diameter
n_dm_uc = 1e4; % number of data
uc2 = @(y) rho_eff_100 * (y / 100) .^ (D_m - 3); % on-demand function...
    % ...for the forward correlation in the mass-mobility domain...
    % ...(rho_eff in [kg/m3] as a function of dm in [nm])
r_uc2 = (dm_lim_uc(2) / dm_lim_uc(1)) ^ (1 / (n_da_uc - 1));
dm_uc = dm_lim_uc(1) * ones(n_da_uc,1) .* r_uc2 .^ (((1 : n_da_uc) - 1)');
rho_eff_uc = uc2(dm_uc);
plt2{end} = plot(dm_uc, rho_eff_uc, 'Color', [0.4940 0.1840 0.5560],...
    'LineStyle', '-.', 'LineWidth', 3);
hold on

% allocate variables to store calculated mobility diameter and effective...
    % ...density
rho_eff = cell(n_dat,1);
dm = cell(n_dat,1);

n_agg = zeros(n_dat,1); % number of aggregates in each timestep saved

rho0 = 1860; % material density

for j = 1 : n_dat
    
    n_agg = length(parsdata(j).npp);

    % initialize effective density array
    rho_eff{j} = zeros(n_agg,1);
    dm{j} = zeros(n_agg,1);
    
    % calculate mobility diameter and effective density
    for jj = 1 : n_agg
        dm{j}(jj) = TRANSP.DIAMOBIL(parsdata(j).dg(jj), parsdata(j).da(jj), fl);
        rho_eff{j}(jj) = rho0 * sum(parsdata(j).pp{jj}(:,2).^3) ./...
            dm{j}(jj).^3;
    end

    % plot temporal data for rho_eff vs. dm
    plt2{j} = scatter(1e9 * dm{j}, rho_eff{j}, ms2(j), mc2(j,:), mt2{j},...
        'LineWidth', 1); 
    
    % make legends and adjust their format
    if j == 1
        legtxt2(j) = strcat('$n_\mathrm{agg}/n_\mathrm{agg_0}$ =',...
            {' '}, num2str(r_n_agg(j), '%.0f'));
    elseif ismember(i, [2,3])
        legtxt2(j) = strcat('$n_\mathrm{agg}/n_\mathrm{agg_0}$ =',...
            {' '}, num2str(r_n_agg(j), '%.1f'));
    else
        legtxt2(j) = strcat('$n_\mathrm{agg}/n_\mathrm{agg_0}$ =',...
            {' '}, num2str(r_n_agg(j), '%.2f'));
    end
end

dx2 = [1e9 * 0.9 * min(cat(1, dm{:})),...
    1e9 * 1.1 * max(cat(1, dm{:}))];
dy2 = [0.9 * min(cat(1, rho_eff{:})), 1.1 * max(cat(1, rho_eff{:}))];

% plot appearance configs
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 11,...
    'TickLength', [0.02 0.02], 'XScale', 'log', 'YScale', 'log')
xlim(dx2)
ylim(dy2)
xlabel('$d_\mathrm{m} [nm]$', 'interpreter', 'latex', 'FontSize', 14)
ylabel('$\rho_\mathrm{eff} [kg/m^3]$', 'interpreter', 'latex', 'FontSize', 14)
box on

% generate legends
lgd2 = legend(cat(1, plt2{:}), legtxt2, 'interpreter', 'latex',...
    'FontSize', 11, 'orientation', 'horizontal', 'NumColumns', 2,...
    'Location', 'southoutside');

% save plots
exportgraphics(f1, 'outputs\dpp-da-temporal.png',...
    'BackgroundColor','none', 'ContentType','vector', 'Resolution', 300)
exportgraphics(f2, 'outputs\rho-dm-temporal.png',...
    'BackgroundColor','none', 'ContentType','vector', 'Resolution', 300)



