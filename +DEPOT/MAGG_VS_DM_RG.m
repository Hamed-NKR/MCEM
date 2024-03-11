function h = MAGG_VS_DM_RG(parsdata, sigmapp_g)
% "PA_VS_NPP" plots the projected area over number of promaries vs. ...
%   ...number of primaries for two sets of hybrid and non-hybrid...
%   ...aggregates and compares them with benchmark values.
% ----------------------------------------------------------------------- %
% 
% Inputs:
%   parsdata: a cell array of structures containing temporal aggregate info
%   sigma_g: geometric standard deviation of agg population pp size
% ----------------------------------------------------------------------- %
%
% Outputs:
%   h: output figure handle
% ----------------------------------------------------------------------- %

% initialize figure 
figure;
h = gcf;
h.Position = [0, 0, 1200, 650];
set(h, 'color', 'white');

% initialize layout
tt = tiledlayout(1,2);
tt.TileSpacing = 'compact';
tt.Padding = 'compact';

% set default for sigma_g
if (~exist('sigma_g', 'var'))  || isempty(sigmapp_g)
    sigmapp_g = [1,1.3]; % aggs are internally monodisperse
end

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

p = cell(n_dat + 2, 2); % initialize the plot cell array
legtxt = cell(n_dat + 2, 1); % placeholder for legends

% set colormap
mc = colormap(hot);
ii = round(1 + (length(mc) - 1) .* (0.05 : 0.7 / (n_dat - 1) : 0.75)');
mc = mc(ii,:);
% mc = flip(mc,1);
mc(6,:) = [236,230,61] / 255;

ms = [25, 25, 25, 45, 27.5, 55]; % Marker sizes
mt = {'o', '^', 'v', 's', 'd', 'p'}; % Marker types

titex = {strcat('(a) $\sigma_{g,pp,ens|_i}$ =', {' '}, num2str(sigmapp_g(1), '%.1f')),...
    strcat('(b) $\sigma_{g,pp,ens|_i}$ =', {' '}, num2str(sigmapp_g(2), '%.1f'))}; % titles for subplots

% reproduce the literature benchmark correlations
% % n0 = logspace(0, 10, 1e4);
% r_dat0 = (1e10 / 1e0)^(1 / (1e4 - 1));
% n0 = 1e0 * ones(1e4,1);
% for i = 2 : 1e4
%     n0(i) = n0(i) * r_dat0^(i-1);
% end
% cor1 = 0.3757 * n0 + 0.4098 * n0.^0.7689;
% 
% cor2 = ((0.94 + 0.03 * sigmapp_g.^4.8) .* (n0.^0.46)).^2; % correlation by Dastanpour & Rogak (2016)

for j = 1 : 2
    nexttile
    
%     % plot correlations
%     p{end - 1, j} = plot(n0, cor1 ./ n0, 'Color', [0.4940 0.1840 0.5560],...
%         'LineStyle', '-.', 'LineWidth', 4);
    hold on
% 
%     p{end, j} = plot(n0, cor2(:,j) ./ n0, 'Color', [0.5 0.5 0.5],...
%         'LineStyle', '--', 'LineWidth', 2.5);
%     
%     if j == 1
%         legtxt{end - 1} = 'Meakin et al. (1989)';
%         legtxt{end} = 'Dastanpour \& Rogak (2016)';
%     end

    % plot temporal aggregate area vs. number data
    for i = 1 : n_dat
        m_agg = zeros(length(parsdata{j}(i).da),1);
        
        for k = 1 : length(parsdata{j}(i).da)
            m_agg(k) = (1860 * pi / 6) * sum(parsdata{j}(i).pp{k}(:,2)).^3;
        end
        
        p{i,j} = scatter(1e9 * 0.75 * parsdata{j}(i).dg, 1e18 * m_agg,...
            ms(i), mc(i,:), mt{i}, 'LineWidth', 1);

        if j == 1
            if i == 1
                legtxt{i} = strcat('$n_{agg}/n_{agg_0}$ =',...
                    {' '}, num2str(t_id(i), '%.0f'));
            else
                legtxt{i} = strcat('$n_{agg}/n_{agg_0}$ =',...
                    {' '}, num2str(t_id(i), '%.2f'));
            end
        end
    end

    % set subplots' properties
    box on
    set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 18,...
        'TickLength', [0.02 0.02], 'XScale', 'log', 'YScale', 'log')
    if j == 2; set(gca, 'yticklabel',[]); end
%     xlim([2, 2e4])
%     ylim([0.4, 0.9])
    
    title(titex{j}, 'FontSize', 22, 'interpreter','latex')
end

% set general plot's properties
xlabel(tt, '$d_m$ [nm]', 'interpreter', 'latex', 'FontSize', 20)
ylabel(tt, '$m_{agg}$ [fg]', 'interpreter', 'latex', 'FontSize', 20)
lgd = legend(cat(2, p{1:end-2,1})', cat(2, legtxt{1:end-2})', 'interpreter', 'latex',...
    'FontSize', 16, 'Orientation', 'horizontal', 'NumColumns', 3);
lgd.Layout.Tile = 'north';    

end

