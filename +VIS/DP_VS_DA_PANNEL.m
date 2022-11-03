function h = DP_VS_DA_PANNEL(parsdata, opts)
% DP_VS_DA_PANNEL makes a pannel of primary particle size with the
%   ...aggregate size data.
% ----------------------------------------------------------------------- %
% 
% Inputs:
%   pardata: a cell array of structures each containing dp and da of...
%       ...aggs over time.
%   opts: plotting options
% ----------------------------------------------------------------------- %
%
% Outputs:
%   h: output figure handle
% ----------------------------------------------------------------------- %

% initialize figure 
figure;
h = gcf;
h.Position = [0, 0, 1200, 900];
set(h, 'color', 'white')
   
% initialize layout
tt = tiledlayout(2,3);
tt.TileSpacing = 'compact';
tt.Padding = 'compact';

% initialize plot & legend placholders
p = cell(7,6);
legtxt = cell(7,1);
legtxt{7} = 'Olfert $\&$ Rogak (2019)';

% reproduce saving criteria for legends
r_kk = (0.02 / 0.5)^(1 / 4);
kk = 0.5 * ones(5,1);
for i = 2 : 5
    kk(i) = kk(i) * r_kk^(i-1);
end
kk = [1; kk];

% set plot title options
if (~exist('opts', 'var')) || (~isfield(opts, 'ttl'))
    opts.ttl = [];
end
titex = opts.ttl;

% set default to the titles
if isempty(titex)
    titex = cell(6,1);
    titex{1} = '(a) Reference';
    titex{2} = '(b) $\overline{d}_{g,agg}$ = 200 nm';
    titex{3} = '(c) $\sigma_{g,agg}$ = 1.6';
    titex{4} = '(d) $\sigma_{g,pp,ens}$ = 1.3';
    titex{5} = '(e) $D_{TEM}$ = 0.5';
    titex{6} = '(f) $VF$ = 20 ppm';
end

% make the colormap and initialize the maker properties
mc = colormap(hot);
ii = round(1 + (length(mc) - 1) .* (0.05 : 0.7 / 5 : 0.75)');
mc = mc(ii,:);
ms = [20, 20, 20, 30, 25, 35];
mt = {'o', '^', 'v', 's', 'd', 'p'};

% generate universal correlation data
dpp_uc = logspace(0, 4, 1000);
da_uc = 100 * (dpp_uc / 17.8).^(1 / 0.35);

for j = 1 : 6
    nexttile
    
    if j ~= 5
        p{7,j} = plot(da_uc, dpp_uc, 'Color', [0.4940 0.1840 0.5560],...
            'LineStyle', '-.', 'LineWidth', 3);
    end
    
    hold on
    
    if contains(titex{j}, 'D_{TEM}')
        legtxt_tem = '$\overline{d}_{pp} = \overline{d}_{pp,100}(d_a/100)^{0.5}$';
        da_tem = 100 * (dpp_uc / 17.8).^(1 / 0.5);
        p_tem = plot(da_tem, dpp_uc, 'Color', mc(1,:), 'LineStyle', ':',...
            'LineWidth', 2.5);
    end    
    
    for i = 1 : 6
        p{i,j} = scatter(1e9 * parsdata{j}(i).da, 1e9 * parsdata{j}(i).dpp_g(:,1),...
            ms(i), mc(i,:), mt{i}, 'LineWidth', 0.1); % plot temporal dp vs. da
        
        % make global legends
        if j == 1
            if i == 1
                legtxt{i} = strcat('$n_{agg}$/$n_{agg_0}$ =',...
                    {' '}, num2str(kk(i), '%.0f'));
            else
                legtxt{i} = strcat('$n_{agg}$/$n_{agg_0}$ =',...
                    {' '}, num2str(kk(i), '%.2f'));
            end
        end
    end
    
    box on
    title(titex{j}, 'FontSize', 24, 'interpreter','latex')        
    xlim([10, 4000])
    ylim([6, 60])
    ax = gca;
    set(ax, 'FontSize', 18, 'TickLength', [0.03 0.03], 'XScale', 'log',...
        'YScale', 'log', 'TickLabelInterpreter','latex')
    ax.XTick = [100,1000];
    if ismember(j, [2,3,5,6])
        set(ax, 'yticklabel',[])
    end
    if ismember(j, [1,2,3])
        set(ax, 'xticklabel',[])
    end    
end

lgd = legend(cat(2, [p{1:6,1}, p_tem, p{7,1}])', cat(2, [legtxt{1:6}, legtxt_tem, legtxt{7}]),...
    'Location', 'northwest', 'Orientation', 'horizontal', 'FontSize', 18,...
    'interpreter', 'latex');
lgd.Layout.Tile = 'north';
lgd.NumColumns = 4;

xlabel(tt, '$d_a$ [nm]', 'FontSize', 20, 'interpreter','latex')
ylabel(tt, '$\overline{d}_{pp}$ [nm]', 'FontSize', 20, 'interpreter', 'latex')

end


