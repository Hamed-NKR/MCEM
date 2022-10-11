function h = NPP_VS_DGDPP_PANNEL(parsdata, opts)
% "PA_VS_NPP" plots the number of primaries vs. gyration diameter over...
%   ...primary particle diameter for two sets of hybrid and non-hybrid...
%   ...aggregates and compares them with benchmark values.
% ----------------------------------------------------------------------- %
% 
% Inputs:
%   parsdata: a cell array of structures containing temporal aggregate info
%   opts: plotting options
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

% initialize plotting options setting
if ~exist('opts', 'var') 
    opts = struct();
end

% field for data fitting and getting fractal properties
if ~isfield(opts, 'fit')
    opts.fit = [];
end
opts_fit = opts.fit;
if isempty(opts_fit)
    opts_fit = 'off'; % default not to fit the data
end

n_dat = 6; % number of time-resolved datasets to be plotted

% data extents
kk_max = 0.5; 
kk_min = 0.02;

r_t = (kk_min / kk_max)^(1 / (n_dat - 2));
t_id = kk_max * ones(n_dat - 1,1);
for i = 2 : (n_dat - 1)
    t_id(i) = t_id(i) * r_t^(i-1);
end
t_id = [1; t_id]; % time identifier of data to be plotted 

if ismember(opts_fit, {'ON', 'On', 'on'})
    p = cell(n_dat + 3,2); % initialize the plot cell array
    legtxt = cell(n_dat + 2,1); % placeholder for legends
else
    p = cell(n_dat + 1,2);
    legtxt = cell(n_dat + 1,1);
end

% set colormap
mc = colormap(hot);
ii = round(1 + (length(mc) - 1) .* (0.05 : 0.7 / (n_dat - 1) : 0.75)');
mc = mc(ii,:);
% mc = flip(mc,1);

ms = [25, 25, 25, 55, 35, 60]; % Marker sizes
mt = {'o', '^', 'v', 's', 'd', 'p'}; % Marker types

titex = {'(a) $\sigma_{g,pp,ens}$ = 1.0',...
    '(b) $\sigma_{g,pp,ens}$ = 1.3'}; % titles for subplots

% plot Sorensen's (2011) recommended benchmark
rd0 = log10(logspace(0, 1e3, 1e4));
n0 = 1.3 * rd0.^1.78;
legtxt{end} = 'Sorensen (2011)';

for j = 1 : 2
    nexttile
    
    p{end,j} = plot(rd0, n0, 'Color', [0.4940 0.1840 0.5560], 'LineStyle', '-.',...
        'LineWidth', 4);
    hold on
    
    for i = 1 : n_dat
        p{i,j} = scatter(parsdata{j}(i).dg ./ parsdata{j}(i).dpp_g(:,1),...
            parsdata{j}(i).npp, ms(i), mc(i,:), mt{i}, 'LineWidth', 0.1); % plot temporal...
                % ...variations of primary number vs gyration-over-primary
                % ...diameter
        
        if i == 1
            legtxt{i} = strcat('$n_{agg}/n_{agg_0}$ =',...
                {' '}, num2str(t_id(i), '%.0f'));
        else
            legtxt{i} = strcat('$n_{agg}/n_{agg_0}$ =',...
                {' '}, num2str(t_id(i), '%.2f'));
        end
        
        if ismember(opts_fit, {'ON', 'On', 'on'}) && (i == n_dat)
            % concatinate the data
            dpp_ens = cat(1, parsdata{j}.dpp_g);
            dpp_ens = dpp_ens(:,1);
            dgg_dpp_ens = cat(1, parsdata{j}.dg) ./ dpp_ens;
            npp_ens = cat(1, parsdata{j}.npp);
            
            fit = fitlm(table(log(dgg_dpp_ens), log(npp_ens)), 'linear',...
                'Weights', sqrt(npp_ens)); % fit a linear regression weighted...
                    % ...by sqrt of number of primaries
            
            df = fit.Coefficients.Estimate(2); % ensemble fractal dimension
            kf = exp(fit.Coefficients.Estimate(1)); % ~ prefactor
            
            % 95% confidence intervals for fractal properties
            ci = coefCI(fit);
            ci_df = ci(2,:);
            ci_kf = exp(ci(1,:));
            
            % 95% ci error bars
            dci_df = max(ci_df) - df;
            dcip_kf = max(ci_kf) - kf;
            dcin_kf = kf - min(ci_kf);
            
            legtxt{end - 1} = 'Ensemble fit';
            
            % generate the fit data
            npp = kf * rd0.^df;
            ci_npp = [ci_kf(1) * rd0.^ci_df(1); ci_kf(2) * rd0.^ci_df(2)];
            
            % plot the main fit and CI bounds
            p{end - 2, j} = plot(rd0, npp, 'Color', [1, 0, 1],...
                'LineStyle', '--', 'LineWidth', 2.5);
            p{end - 1, j} = plot(rd0, ci_npp(1,:), rd0, ci_npp(2,:),...
                'Color', [1, 0, 1], 'LineStyle', ':', 'LineWidth', 2);
        end
    end
    
    % set subplots' properties
    box on
    set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 18,...
        'TickLength', [0.02 0.02], 'XScale', 'log', 'YScale', 'log')
    if j == 2; set(gca, 'yticklabel',[]); end
    xlim([2, 3e2])
    ylim([1, 2e4])
    
    title(titex{j}, 'FontSize', 22, 'interpreter','latex')
    
    % print fractal properties
    if j == 1
        text(6, 4.5e3, '$d_{f_0} = 1.78, k_{f_0} = 1.3$', 'interpreter', 'latex',...
            'FontSize', 16, 'Color', [0.4940 0.1840 0.5560])
        annotation('arrow', [0.37,0.32], [0.71,0.71], 'Color', [0.4940 0.1840 0.5560])
    end
    text(8, 6, strcat('$d_{f_{ens}}$ =', {' '}, num2str(df, '%.2f'), {' '},...
        '$\pm$', {' '}, num2str(dci_df, '%.2f'), {','}, string(newline),...
        '  $k_{f_{ens}}$ =', {' '}, num2str(kf, '%.2f'), {' '},...
        {'$\pm$'}, {' '}, num2str(max(dcip_kf, dcin_kf), '%.2f')),...
        'interpreter', 'latex', 'FontSize', 16, 'Color', [1, 0, 1])
    if j == 1
        annotation('arrow', [0.122,0.172], [0.27,0.27], 'Color', [1, 0, 1])
    else
        annotation('arrow', [0.598,0.648], [0.27,0.27], 'Color', [1, 0, 1])
    end
    
%     text(7, 6, strcat('$d_{f_{ens}}$ =', {' '}, num2str(df, '%.2f'),...
%         ', 95$\%$ CI = [+', num2str(dci_df, '%.2f'), {',-'}, num2str(dci_df, '%.2f'),...
%         {'],'}, string(newline), '  $k_{f_{ens}}$ =', {' '}, num2str(kf, '%.2f'),...
%         {', 95$\%$ CI = [+'}, num2str(dcip_kf, '%.2f'), {',-'}, num2str(dcin_kf, '%.2f'),...
%         {']'}), 'interpreter', 'latex', 'FontSize', 16,...
%         'Color', [1, 0, 1])
end

% set general plot's properties
xlabel(tt, '$d_g/d_{g,pp}$ [-]', 'interpreter', 'latex', 'FontSize', 20)
ylabel(tt, '$n_{pp}$ [-]', 'interpreter', 'latex', 'FontSize', 20)
if ismember(opts_fit, {'ON', 'On', 'on'})
    lgd = legend(cat(2, [p{1 : n_dat, 1}, p{n_dat + 1, 1}, p{n_dat + 3, 1}])',...
        cat(2, [legtxt{1 : n_dat}, legtxt{n_dat + 1}, legtxt{n_dat + 2}])',...
        'interpreter', 'latex', 'FontSize', 18, 'Orientation', 'horizontal',...
        'NumColumns', 4);
else
    lgd = legend(cat(2, p{1 : n_dat + 1, 1})', cat(2, legtxt{1 : n_dat + 1})',...
        'interpreter', 'latex', 'FontSize', 18, 'Orientation', 'horizontal',...
        'NumColumns', 4);
end
lgd.Layout.Tile = 'north';

end

% 
% $d_f = 1.78, k_f = 1.3$'