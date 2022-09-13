function h = PA_VS_NPP(parsdata)
% "PA_VS_NPP" plots the projected area over number of promaries vs. ...
%   ...number of primaries for two sets of hybrid and non-hybrid...
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
h.Position = [0, 0, 800, 800];
set(h, 'color', 'white');

n_dat = 5; % number of hybrid datasets to be plotted

% data extents
kk_max = 0.5; 
kk_min = 0.02;

r_t = (kk_min / kk_max)^(1 / (n_dat - 1));
t_id = kk_max * ones(n_dat,1);
for i = 2 : n_dat
    t_id(i) = t_id(i) * r_t^(i-1);
end
t_id = [1; t_id]; % time identifier of data to be plotted 

p = cell(n_dat + 2,1); % initialize the plot cell array
legtxt = cell(n_dat + 2,1); % placeholder for legends

% set colormap
mc = colormap(hot);
ii = round(1 + (length(mc) - 1) .* (0.05 : 0.8 / n_dat : 0.85)');
mc = mc(ii,:);
% mc = flip(mc,1);

ms = [25, 25, 25, 55, 35, 60]; % Marker sizes
mt = {'o', '^', 'v', 's', 'd', 'p'}; % Marker types

% plot the Meakin et al. (1989) recommended correlation
% n0 = log(logspace(exp(1), exp(10), 1e4));
r_dat0 = (1e10 / 1e0)^(1 / (1e4 - 1));
n0 = 1e0 * ones(1e4,1);
for i = 2 : 1e4
    n0(i) = n0(i) * r_dat0^(i-1);
end
a0 = 0.3757 * n0 + 0.4098 * n0.^0.7689;
p{1} = plot(n0, a0 ./ n0, 'Color', [0.4940 0.1840 0.5560], 'LineStyle', '-.',...
    'LineWidth', 2.5);
legtxt{1} = 'Meakin et al. (1989)';
hold on

% plot temporal aggregate area vs. number data
for i = 1 : (n_dat + 1)
    p{i + 1} = scatter(parsdata(i).npp, ((parsdata(i).da ./...
        parsdata(i).dpp_g(:,1)).^2 * pi / 4) ./ parsdata(i).npp,...
        ms(i), mc(i,:), mt{i}, 'LineWidth', 0.1);
    
    legtxt{i + 1} = strcat('\itn_{agg}\rm/\itn_{agg_0}\rm =',...
        {' '}, num2str(t_id(i), '%.2f'));
end

% set axes
box on
set(gca, 'FontName', 'Calibri Light', 'FontSize', 12, 'TickLength', [0.005 0.005])
xlabel('{\itn_{pp}} [-]', 'FontName', 'Calibri', 'FontSize', 16)
xlim([2, 2e4])
set(gca, 'XScale', 'log')
ylabel('{\itA_{agg}}/{\itn_{pp}} [-]', 'FontName', 'Calibri', 'FontSize', 16)
ylim([-inf, 1])
set(gca, 'YScale', 'log')
legend([p{1}, cat(1, p{2:end})'], cat(2, legtxt{:}), 'FontName', 'Calibri Light',...
    'FontSize', 12, 'Location', 'northeast')

end

