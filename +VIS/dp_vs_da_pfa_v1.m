function [h, D_pol] = dp_vs_da_pfa_v1(parsdata, params, opts)
% "dp_vs_da_pfa_v1" plots a group of primary particle size vs aggregate...
%   ...size plots within each the evolution of aggregates from a...
%   ...totally monodisperse toward a hybrid structure with nonuniform...
%   ...primaries is presented. The plots vary depending on the...
%   ...initialization parameters.
% ----------------------------------------------------------------------- %
% 
% Inputs:
%   pardata: a cell array of structures each containing dp, da and stds of
%       ...aggs over time staring from the monodispersity moment.
%   params: a structure containing parameters to be plotted.
%   opts: plotting options
% ----------------------------------------------------------------------- %
%
% Outputs:
%   D_pol: The scaling exponent series of dp vs. da data
%   h: output figure handle
% ----------------------------------------------------------------------- %

% initialize figure 
figure;
h = gcf;
h.Position = [0, 0, 1800, 600];
set(h, 'color', 'white');

% initialize layout
tt = tiledlayout(1,3);
tt.TileSpacing = 'compact';
tt.Padding = 'compact';

% initialize plot & legend placholders
p = cell(13,3);
D_pol = zeros(6,3);
legtxt = cell(7,4);
legtxt{7,4} = 'Olfert & Rogak (2019)';

% reproduce saving criteria for legends
r_kk = (0.02 / 0.5)^(1 / 4);
kk = 0.5 * ones(5,1);
for i = 2 : 5
    kk(i) = kk(i) * r_kk^(i-1);
end
kk = [1; kk];

% the independent variable to be used for polydispersity exponenet fitting
da_fit = logspace(0, 4, 1000);

% initialize maker properties
mc = colormap(turbo); % make colormap to monitor hybridization
ii = round(1 + (length(mc) - 1) .* (0.05 : 0.9 / 5 : 0.95)'); % define the sample points over colormap range
mc = mc(ii,:); % sample across the colormap for different lifestages
mc = flip(mc,1); % reverse the colormap direction (if needed)
ms = [25, 25, 25, 35, 35, 50]; % Marker sizes
mt = {'o', '^', 'v', 's', 'd', 'p'}; % Marker type

% generate universal correlation data
dpp_uc = log(logspace(1, 100, 1000));
da_uc = 100 * (dpp_uc / 17.8).^(1 / 0.35);

for j = 1 : 3
    nexttile
    
    for i = 1 : 6
        p{2*i-1,j} = scatter(1e9 * parsdata{j}(i).da, 1e9 * parsdata{j}(i).dpp_g(:,1),...
            ms(i), mc(i,:), mt{i}, 'LineWidth', 0.1); % plot temporal dp vs. da
        hold on
        
        % make global legends
        if j == 1
            if i == 1
                legtxt{i,4} = strcat('\itn_{agg}\rm/\itn_{agg_0}\rm =',...
                    {' '}, num2str(kk(i), '%.0f'));
            else
                legtxt{i,4} = strcat('\itn_{agg}\rm/\itn_{agg_0}\rm =',...
                    {' '}, num2str(kk(i), '%.2f'));
            end
        end
        
        fit_pol = fitlm(table(log(parsdata{j}(i).da), log(parsdata{j}(i).dpp_g(:,1))), 'linear'); % find a fit by linear regression
        D_pol(i,j) = fit_pol.Coefficients.Estimate(2); % get polydispersity exponent from the fit
        dpp_fit = 1e9 * exp(D_pol(i,j) * log(1e-9 * da_fit) + fit_pol.Coefficients.Estimate(1)); % generate the fit data
        p{2*i,j} = plot(da_fit, dpp_fit, 'Color', mc(i,:),...
            'LineStyle', '--', 'LineWidth', 1); % plot the fit
        
        legtxt{i,j} = strcat('D_{pol} =', {' '}, num2str(D_pol(i,j), '%.2f')); % fit legends
    end
    
    legtxt{7,j} = 'D_{pol_{uc}} = 0.35'; % legend of universal correlation
    
    p{13,j} = plot(da_uc, dpp_uc, 'Color', [0.5 0.5 0.5],...
        'LineStyle', '-.', 'LineWidth', 2.5); % plot universal correlation
    
    box on
    set(gca, 'FontName', 'SansSerif', 'FontSize', 12, 'TickLength', [0.02 0.02])
    xlabel('{\it d_a} [nm]', 'FontName', 'SansSerif', 'FontSize', 14)
    xlim([5, 5000])
    set(gca, 'XScale', 'log')
    if j == 1
        ylabel('{\it d_{pp}} [nm]', 'FontName', 'SansSerif', 'FontSize', 14)
    end
    ylim([5, 100])
    set(gca, 'YScale', 'log')
    ttltxt = strcat('{\it', {' '}, params.name, '}', '=', {' '},...
        num2str(params.value{j}, '%.0f'), {' '}, '[', params.unit,...
        ']');
    title(ttltxt, 'FontName', 'SansSerif', 'FontWeight', 'bold', 'FontSize', 16)
    
    if j == 3
        p_dummy = plot(da_uc, dpp_uc, 'Color', [0.5 0.5 0.5],...
            'LineStyle', '-.', 'LineWidth', 2.5);
        legend(cat(1, p{1:2:13,3})', cat(2, legtxt{:,4}),...
            'FontName', 'SansSerif', 'FontSize', 12,...
            'Location', 'eastoutside');
        ax = axes('position', get(gca,'position'), 'visible', 'off');
        legend(ax, [cat(1, p{2:2:12,j})', p_dummy], cat(2, legtxt{:,j}),...
            'Location', 'northwest', 'FontName', 'SansSerif', 'FontSize', 10);        
    else
        legend([cat(1, p{2:2:12,j})', p{13,j}], cat(2, legtxt{:,j}),...
            'Location', 'northwest', 'FontName', 'SansSerif', 'FontSize', 10);
    end    
end

