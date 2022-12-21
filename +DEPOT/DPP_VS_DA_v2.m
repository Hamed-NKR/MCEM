function h = DPP_VS_DA_v2(parsdata)
% DPP_VS_DA_v2 makes a temporal scatter graph of primary particle size vs...
%   ...aggregate size variation data.
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
h.Position = [0, 0, 700, 650];
set(h, 'color', 'white')
   
% initialize plot & legend placholders
p = cell(7,1);
legtxt = cell(7,1);
legtxt{7} = 'Olfert $\&$ Rogak (2019)';

% reproduce saving criteria for legends
r_kk = (0.02 / 0.5)^(1 / 4);
kk = 0.5 * ones(5,1);
for i = 2 : 5
    kk(i) = kk(i) * r_kk^(i-1);
end
kk = [1; kk];

% make the colormap and initialize the maker properties
mc = colormap(hot);
ii = round(1 + (length(mc) - 1) .* (0.05 : 0.7 / 5 : 0.75)');
mc = mc(ii,:);
ms = [20, 20, 20, 30, 25, 35];
mt = {'o', '^', 'v', 's', 'd', 'p'};

% generate universal correlation data
dpp_uc = logspace(0, 4, 1000);
da_uc = 100 * (dpp_uc / 17.8).^(1 / 0.35);

p{7,1} = plot(da_uc, dpp_uc, 'Color', [0.4940 0.1840 0.5560],...
    'LineStyle', '-.', 'LineWidth', 3);
hold on

for i = 1 : 6
    p{i,1} = scatter(1e9 * parsdata(i).da, 1e9 * parsdata(i).dpp_g(:,1),...
        ms(i), mc(i,:), mt{i}, 'LineWidth', 0.1); % plot temporal dp vs. da
    
    % make global legends
    if i == 1
        legtxt(i) = strcat('$n_{agg}$/$n_{agg_0}$ =',...
            {' '}, num2str(kk(i), '%.0f'));
    else
        legtxt(i) = strcat('$n_{agg}$/$n_{agg_0}$ =',...
            {' '}, num2str(kk(i), '%.2f'));
    end
end

box on
xlim([35, 2000])
ylim([10, 40])
ax = gca;
set(ax, 'FontSize', 18, 'TickLength', [0.02 0.02], 'XScale', 'log',...
    'YScale', 'log', 'TickLabelInterpreter','latex')
ax.XTick = [100, 1000];
ax.YTick = [10,20, 30, 40, 50];

legend(cat(1, p{:})', cat(2,legtxt(:)), 'Location', 'northwest',...
    'FontSize', 16, 'interpreter', 'latex');

xlabel('$d_a$ [nm]', 'FontSize', 20, 'interpreter','latex')
ylabel('$\overline{d}_{pp}$ [nm]', 'FontSize', 20, 'interpreter', 'latex')

end


