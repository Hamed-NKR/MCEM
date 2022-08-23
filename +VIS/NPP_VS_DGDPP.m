function h = NPP_VS_DGDPP(parsdata)
% "PA_VS_NPP" plots the number of primaries vs. gyration diameter over...
%   ...primary particle diameter for two sets of hybrid and non-hybrid...
%   ...aggregates and compares them with benchmark values.
% ----------------------------------------------------------------------- %
% 
% Inputs:
%   parsdata: a cell array of structures containing temporal aggregate info
% ----------------------------------------------------------------------- %
%
% Outputs:
%   h: output figure handle
% ----------------------------------------------------------------------- %

% initialize figure 
figure;
h = gcf;
h.Position = [0, 0, 1800, 600];
set(h, 'color', 'white');

r_t = (0.02 / 0.5)^(1 / 4);
t_id = 0.5 * ones(5,1);
for i = 2 : 5
    t_id(i) = t_id(i) * r_t^(i-1);
end
t_id = [1; t_id]; % time identifier of data to be plotted 

p = cell(7,1); % initialize the plot cell array
legtxt = cell(7,1); % placeholder for legends

% set colormap
mc = colormap(hot);
ii = round(1 + (length(mc) - 1) .* (0.05 : 0.7 / 5 : 0.75)');
mc = mc(ii,:);
% mc = flip(mc,1);

ms = [25, 25, 25, 35, 35, 50]; % Marker sizes
mt = {'o', '^', 'v', 's', 'd', 'p'}; % Marker types

% plot the Meakin et al. (1989) recommended correlation
n0 = log(logspace(1, 1e5, 1e6));
rd0 = exp((log(n0) - 1.3) / 1.78);
p{1} = plot(rd0, n0, 'Color', [0.5 0.5 0.5], 'LineStyle', '-.',...
    'LineWidth', 2.5);
legtxt{1} = 'Sorensen et al., 2011';

% plot temporal aggregate area vs. number data
for i = 1 : 6
    p{i + 1} = scatter(parsdata(i).dg ./ parsdata(i).dpp_g,...
        parsdata(i).npp, ms(i), mc(i,:), mt{i}, 'LineWidth', 0.1);
    
    legtxt{i + 1} = strcat('\itn_{agg}\rm/\itn_{agg_0}\rm =',...
        {' '}, num2str(t_id(i), '%.2f'));
end

% set axes
box on
set(gca, 'FontName', 'Calibri Light', 'FontSize', 12, 'TickLength', [0.005 0.005])
xlabel('{\itn_{pp}} [-]', 'FontName', 'Calibri', 'FontSize', 16)
set(gca, 'XScale', 'log')
ylabel('{\itd_g}/{\itd_{g,pp}} [-]', 'FontName', 'Calibri', 'FontSize', 16)
set(gca, 'YScale', 'log')
legend(cat(1, p(:)), cat(2, legtxt{:}), 'FontName', 'Calibri Light', 'FontSize', 12,...
    'Location', 'northwest')

end

