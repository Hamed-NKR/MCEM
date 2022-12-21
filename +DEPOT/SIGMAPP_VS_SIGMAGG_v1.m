function h = SIGMAPP_VS_SIGMAGG_v1(parsdata, params, opts)
% "SIGMAPP_VS_SIGMAGG" plots boxcharts from temproal evolution of...
%   ...standard deviation of primary particle size within vs. bewteen...
%   ...the aggregates.
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
%   h: output figure handle
% ----------------------------------------------------------------------- %

% initialize figure 
figure;
h = gcf;
h.Position = [0, 0, 1800, 600];
set(h, 'color', 'white');

% set data fitting option
if (~exist('opts', 'var')) || (~isfield(opts, 'scat'))
    opts.scat = [];
end
opts_scat = opts.scat;
if isempty(opts_scat)
    opts_scat = 'off'; % default not to plot the scatter points
end

tbl = []; % initialize the datasheet to be plotted

n_dat = zeros(6,3); % number of data points in each boxplot

r_t = (0.02 / 0.5)^(1 / 4);
t_id = 0.5 * ones(5,1);
for i = 2 : 5
    t_id(i) = t_id(i) * r_t^(i-1);
end
t_id = [1; t_id]; % time identifier of data to be plotted 

legtxt = cell(6,1); % placeholder for legends

% set box and point colors
mc = colormap(hot); % choose colormap
ii = round(1 + (length(mc) - 1) .* (0.05 : 0.7 / 5 : 0.75)'); % define the sample points over colormap range
mc = mc(ii,:); % sample across the colormap for different lifestages
% mc = flip(mc,1); % reverse the colormap direction (if needed)

% extract data from inputs
for j = 1 : 3
    for i = 1 : 6
        n_dat(i,j) = size(parsdata{j}(i).dpp_g, 1);
        tbl_temp = [repmat(t_id(7 - i), n_dat(i,j), 1),...
            repmat(params.var(j), n_dat(i,j), 1), parsdata{j}(i).dpp_g(:,2)];
        tbl = vertcat(tbl, tbl_temp);
        
        if j == 1
            legtxt{i} = strcat(' \itn_{agg}\rm/\itn_{agg_0}\rm =',...
                {' '}, num2str(t_id(i), '%.2f'));
        end
    end
end

x_bc = categorical(tbl(:,2)); % generate the plot x axis

% make the boxplot
bc = boxchart(x_bc, tbl(:,3), 'GroupByColor', tbl(:,1));
for i = 1 : 6
    bc(i).BoxFaceColor = mc(i,:);
    bc(i).MarkerColor = mc(i,:);
end
hold on

dx = 1/6; % gap to consider between the boxes to locate the scatter data

% put scatter data on top of boxplots
if ismember(opts_scat, {'ON', 'On', 'on'})
    for j = 1 : 3
        for i = 1 : 6
            x_scat = repmat(j + ((2 * i - 1) / 2 - 3) * dx, n_dat(i,j), 1);
            scatter(x_scat, parsdata{j}(i).dpp_g(:,2), [], mc(i,:),...
                'MarkerEdgeAlpha', 0.3)
        end
    end
end

box on
set(gca, 'FontName', 'Calibri Light', 'FontSize', 12, 'TickLength', [0.005 0.005])
xlabel('{\it {\sigma}_{g,pp,ens}} [-]', 'FontName', 'Calibri', 'FontSize', 16)
ylabel('{\it {\sigma}_{g,pp,agg}} [-]', 'FontName', 'Calibri', 'FontSize', 16)
ylim([1, ceil(max(ylim) / 0.1) * 0.1])
legend(bc, cat(2, legtxt{:}), 'FontName', 'Calibri Light', 'FontSize', 12,...
    'Location', 'northeastoutside')

end

