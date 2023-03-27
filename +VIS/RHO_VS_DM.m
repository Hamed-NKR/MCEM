function h = RHO_VS_DM(parsdata)
% "RHO_VS_DM" plots the effective density vs. mobility diameter...
%   ...for a certain polydispersity levels.
% ----------------------------------------------------------------------- %
% 
% Input:
%   parsdata: a cell array of structures containing temporal aggregate info
% ----------------------------------------------------------------------- %
%
% Output:
%   h: output figure handle
% ----------------------------------------------------------------------- %

% initialize figure 
figure;
h = gcf;
h.Position = [0, 0, 700, 700];
set(h, 'color', 'white')

% data saving moments
n_dat = 6;
kk_max = 0.5; 
kk_min = 0.02;
r_t = (kk_min / kk_max)^(1 / (n_dat - 2));
t_id = kk_max * ones(n_dat - 1,1);
for i = 2 : (n_dat - 1)
    t_id(i) = t_id(i) * r_t^(i-1);
end
t_id = [1; t_id];

% setting properties for guidelines of mass-mobility exponent
dm_gl = [1.75, 1.9, 2.05, 2.2, 2.35];
km_gl = [800, 700, 650, 600, 550];
n_gl = length(dm_gl); % number of guidelines for effective density
x_gl = [1.25, 1.25, 1.25, 1.25, 1.25]; % x location of guideline label
y_gl = [275000, 120000, 57000, 27000, 12500]; % y location ~
theta_gl = [42, 39, 35, 31, 26]; % orientation of ~

p = cell(n_dat + n_gl + 1, 1); % initialize the plot cell array
legtxt = cell(n_dat + n_gl + 1, 1); % placeholder for legends

% set colormap
mc = colormap(hot);
ii = round(1 + (length(mc) - 1) .* (0.05 : 0.7 / (n_dat - 1) : 0.75)');
mc = mc(ii,:);
% mc = flip(mc,1);
mc(6,:) = [236,230,61] / 255;

ms = [25, 25, 25, 45, 27.5, 55]; % Marker sizes
mt = {'o', '^', 'v', 's', 'd', 'p'}; % Marker types

% reproduce the literature benchmark correlations
% n0 = logspace(0, 10, 1e4);
r0 = (1e4 / 1e0)^(1 / (1e4 - 1));
dm0 = 1e0 * ones(1e4,1);
for i = 2 : 1e4
    dm0(i) = dm0(i) * r0^(i-1);
end
rho0 = 510 * (dm0 / 100).^(-0.52);

% plot universal correlation
p{end - 5} = plot(dm0, rho0, 'Color', [0.4940 0.1840 0.5560],...
    'LineStyle', '-.', 'LineWidth', 4);

hold on

% label unversal correlation
legtxt{end - 5} = 'Olfert and Rogak (2019)';
t0 = text(1.5, 2500, strcat('\boldmath$d_m =$', {' \bf'}, num2str(2.48)),...
    'interpreter', 'latex', 'FontSize', 14, 'Color', [0.4940 0.1840 0.5560]);
set(t0, 'Rotation', -19);

% plot guidelines and label them
for ii = 1 : n_gl
    rho_gl = km_gl(ii) * (dm0 / 100).^(dm_gl(ii) - 3);
    p{end - ii + 1} = plot(dm0, rho_gl, 'Color', [0.1 0.1 0.1],...
        'LineStyle', '--', 'LineWidth', 1.5);
    
    t_gl = text(x_gl(ii), y_gl(ii), strcat('$d_m =$', {' '}, num2str(dm_gl(ii), '%.2f')),...
        'interpreter', 'latex', 'FontSize', 12, 'Color', [0.1 0.1 0.1]);
    set(t_gl, 'Rotation', -theta_gl(ii));    
end

% plot temporal effective density vs. mobility diameter data
for i = 1 : n_dat
    n_agg_i = length(parsdata(i).da);
    
    rho_eff_i = zeros(n_agg_i,1);
    
    dm_i = TRANSP.DIAMOBIL(parsdata(i).dg, parsdata(i).da);
    
    for k = 1 : n_agg_i
        rho_eff_i(k) = 1860 * sum(parsdata(i).pp{k}(:,2).^3) ./ dm_i(k).^3;
    end
    
    p{i} = scatter(1e9 * dm_i, rho_eff_i, ms(i), mc(i,:), mt{i}, 'LineWidth', 1);
    
    if i == 1
        legtxt{i} = strcat('$n_{agg}/n_{agg_0}$ =',...
            {' '}, num2str(t_id(i), '%.0f'));
    else
        legtxt{i} = strcat('$n_{agg}/n_{agg_0}$ =',...
            {' '}, num2str(t_id(i), '%.2f'));
    end
end

% set subplots' properties
box on
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 18,...
    'TickLength', [0.02 0.02], 'XScale', 'log', 'YScale', 'log')
% xlim([10, 1e4])
% ylim([1, 1e3])

% set general plot's properties
xlabel('$d_m$ [nm]', 'interpreter', 'latex', 'FontSize', 20)
ylabel('$\rho_{eff} [kg/m^{3}]$', 'interpreter', 'latex', 'FontSize', 20)
legend(cat(2, p{1 : end - 5})', cat(2, legtxt{1 : end - 5})', 'interpreter', 'latex',...
    'FontSize', 16, 'location', 'northeast')

end

