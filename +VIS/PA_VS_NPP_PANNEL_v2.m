function h = PA_VS_NPP_PANNEL_v2(parsdata, sigmapp_g)
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

% % set colormap
% mc = colormap(hot);
% ii = round(1 + (length(mc) - 1) .* (0.05 : 0.7 / (n_dat - 1) : 0.75)');
% mc = mc(ii,:);
% % mc = flip(mc,1);
% mc(6,:) = [236,230,61] / 255;

ms = [25, 25, 25, 45, 27.5, 55]; % Marker sizes
mt = {'o', '^', 'v', 's', 'd', 'p'}; % Marker types

titex = {strcat('(a) $\sigma_{g,pp,ens|_i}$ =', {' '}, num2str(sigmapp_g(1), '%.1f')),...
    strcat('(b) $\sigma_{g,pp,ens|_i}$ =', {' '}, num2str(sigmapp_g(2), '%.1f'))}; % titles for subplots

% reproduce the literature benchmark correlations
% n0 = logspace(0, 10, 1e4);
r_dat0 = (1e10 / 1e0)^(1 / (1e4 - 1));
n0 = 1e0 * ones(1e4,1);
for i = 2 : 1e4
    n0(i) = n0(i) * r_dat0^(i-1);
end
cor1 = 0.3757 * n0 + 0.4098 * n0.^0.7689;

% set guideline properties
sigma_gl = 1 : 0.1 : 1.3;
n_gl = length(sigma_gl);
cor2 = ((0.94 + 0.03 * sigma_gl.^4.8) .* (n0.^0.46)).^2; % guidlines from...
    % ...correlation by Dastanpour & Rogak (2016)
x_gl = [3, 9.25, 13.5, 25]; % x location of guideline labels
y_gl = [0.8225, 0.8375, 0.855, 0.87]; % y location ~
theta_gl = 43; % orientation of ~

for j = 1 : 2
    nexttile
    
    % plot correlations
    p{end - 1, j} = plot(n0, cor1 ./ n0, 'Color', [0.4940 0.1840 0.5560],...
        'LineStyle', '-.', 'LineWidth', 4);
    hold on
    
    for ii = 1 : n_gl
        if ii == 1
            p{end, j} = plot(n0, cor2(:,ii) ./ n0, 'Color', [0.5 0.5 0.5],...
                'LineStyle', '--', 'LineWidth', 2.5);
            t_gl = text(x_gl(ii), y_gl(ii), strcat('$\sigma_{g,pp} =$', {' '},...
                 num2str(sigma_gl(ii), '%.1f'), ', from', string(newline), ' Sorensen (2011)'),...
                 'interpreter', 'latex', 'FontSize', 12, 'Color', [0.1 0.1 0.1]);
            set(t_gl, 'Rotation', -theta_gl);
        else
            plot(n0, cor2(:,ii) ./ n0, 'Color', [0.5 0.5 0.5],...
                'LineStyle', '--', 'LineWidth', 2.5);
            t_gl = text(x_gl(ii), y_gl(ii), strcat('$\sigma_{g,pp} =$', {' '},...
                 num2str(sigma_gl(ii), '%.1f')),'interpreter', 'latex',...
                 'FontSize', 12, 'Color', [0.1 0.1 0.1]);
            set(t_gl, 'Rotation', -theta_gl);
        end
    end
    
    if j == 1
        legtxt{end - 1} = 'Meakin et al. (1989)';
        legtxt{end} = 'Dastanpour \& Rogak (2016)';
    end

    % plot temporal aggregate area vs. number data
    for i = 1 : n_dat
        dpp0 = zeros(length(parsdata{j}(i).npp), 1);
        for k = 1 : length(dpp0)
            pp0 = parsdata{j}(i).pp{k};
            dpp0(k) = sqrt(sum(pp0(:,2)).^2) / parsdata{j}(i).npp(k);
        end
                
        p{i,j} = scatter(parsdata{j}(i).npp, ((parsdata{j}(i).da ./...
            dpp0).^2 * pi / 4) ./ parsdata{j}(i).npp,...
            ms(i), parsdata{j}(i).dpp_g(:,2), mt{i}, 'LineWidth', 1);
        
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
    xlim([2, 2e4])
    ylim([0.4, 0.9])
    
    title(titex{j}, 'FontSize', 22, 'interpreter','latex')
end

cb = colorbar;
colormap turbo
cb.Label.String = '$\sigma_{g,pp,agg}$ [-]';
cb.Label.Interpreter  = 'latex';
cb.Label.Rotation = 360;
cb.TickLabelInterpreter  = 'latex';
cb.FontSize = 16;
cb.Label.FontSize = 20;
cb.Limits = [1.0 1.7];
cb.Ticks = 1.0:0.1:1.7;
cb.TickLength = 0.02;
cb.Layout.Tile = 'east';
cbpos = get(cb, 'Position');
cb.Label.Position = [cbpos(1) - 0.75 , cbpos(2) + 1.705];


% dummy plots for the legend appearance purposes
pnot = cell(n_dat,1);
nppnot = 1e-3 * (1 : n_dat);
Abarnot = 1e-3 * (1 : n_dat);
for i = 1 : n_dat
    pnot{i} = scatter(nppnot(i), Abarnot(i), ms(i), [0.1 0.1 0.1], mt{i});
end

% set general plot's properties
xlabel(tt, '$n_{pp}$ [-]', 'interpreter', 'latex', 'FontSize', 20)
ylabel(tt, '$\hat{\overline{A}}_{agg}/n_{pp}$ [-]', 'interpreter', 'latex', 'FontSize', 20)
lgd = legend(cat(2, [pnot{:}, p{n_dat + 1:end,1}])', cat(2, legtxt{:})', 'interpreter', 'latex',...
    'FontSize', 16, 'Orientation', 'horizontal', 'NumColumns', 4);
lgd.Layout.Tile = 'north';    

end

