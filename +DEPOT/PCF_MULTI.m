function h = PCF_MULTI(g, r, n_hyb_0, sigma_pp_0)
% PCF_MULTI gets multiple pair-correlation function (PCF) distributions...
%   ...from different aggregates/scenarios and visualizes them for...
%   ...comparison in a single PCF vs. radial location plot. 
% ----------------------------------------------------------------------- %
% 
% Inputs:
%   g: A cell array containing data of radial vector of averaged (over...
%       ...primary particles) PCF within and nearby different aggregates.
%   r: The raidal logarithmic points of calculation for PCF starting...
%       ...from inside a central primary particle. This can also be...
%       ...different for each aggregates.
%   n_hyb_0: A vector containing number of sub-agggregates for individual
%       ...aggregates
%   sigma_pp_0: The geometric standard deviation of primary particle size
%       ...withing aggregates
% ----------------------------------------------------------------------- %
% 
% Outputs:
%   h: PCF plot figure handle
% ----------------------------------------------------------------------- %

% initialize figure 
figure;
h = gcf;
h.Position = [0, 0, 800, 500];
set(h, 'color', 'white');

% a placeholder for unique values of number of subaggregates
n_hyb_00 = unique(n_hyb_0);
% count_hyb = length(n_hyb_00);

n_agg = length(g); % number of aggregates to be plotted

plt = cell(n_agg, 1); % initialize plot variable
lbl = cell(n_agg, 1); % initialize placeholder for legends

lin_c = [252,141,98;
    102,194,165;
    141,160,203] / 255; % line color repository
lin_w = [1.5, 1.5, 1.5]; % line width repo.
lin_t = {'-', '--', ':', '-.', '--'}; % line style repo.

for i = 1 : n_agg
    ii = find(n_hyb_00 == n_hyb_0(i)); % label aggregates with their number of sub-aggregates
    iii = i - find(n_hyb_0 == n_hyb_00(ii),1) + 1; % identify the order of datapoints
    
    plt{i} = plot(r{i}, g{i}, 'Color', lin_c(ii,:), 'LineStyle', lin_t{iii},...
        'LineWidth', lin_w(ii)); % plot pair-correlation function
    if iii > 4; plt{i}.Marker = '+'; plt{i}.MarkerSize = 5; set(plt{i}, 'linewidth', 0.5); end
    
    % make the legend
    lbl{i} = strcat({'$n_\mathrm{hyb}$ ='}, {' '}, num2str(n_hyb_0(i), '%d'),...
        {','}, ' $\sigma_\mathrm{g,pp}$ =', {' '}, num2str(sigma_pp_0(i), '%.2f'));    
    
    hold on
end

box on
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 18,...
    'TickLength', [0.02 0.02], 'XScale', 'log', 'YScale', 'log')
xlabel('$\overline{r}$ (-)', 'interpreter', 'latex', 'FontSize', 20)
ylabel('$\overline{g}$($\overline{r}$) (-)', 'interpreter', 'latex', 'FontSize', 20)
xlim([0.05 120])
ylim([5e-5 1])
legend(cat(2, plt{:})', cat(2, lbl{:})', 'interpreter', 'latex', 'FontSize', 16, 'Location', 'eastoutside',...
    'NumColumns', 1);

end

