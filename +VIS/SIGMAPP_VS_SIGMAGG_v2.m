function h = SIGMAPP_VS_SIGMAGG_v2(parsdata, opts)
% "SIGMAPP_VS_SIGMAGG" plots boxcharts from temproal evolution of...
%   ...standard deviation of primary particle size within vs. bewteen...
%   ...the aggregates.
% ----------------------------------------------------------------------- %
%
% Inputs:
%   pardata: a cell array of structures each containing dp, da and stds of...
%       ...aggs over time staring from the monodispersity moment.
%   opts: plotting options
% ----------------------------------------------------------------------- %
%
% Outputs:
%   h: output figure handle
% ----------------------------------------------------------------------- %

% initialize figure 
figure;
h = gcf;
h.Position = [0, 0, 1200, 600];
set(h, 'color', 'white');

% initialize options
if (~exist('opts', 'var'))
    opts = struct();
end

% set data fitting option
if (~isfield(opts, 'scat'))
    opts.scat = [];
end
opts_scat = opts.scat;
if isempty(opts_scat)
    opts_scat = 'off'; % default not to plot the scatter points
end

% set title option
if (~isfield(opts, 'ttl'))
    opts.ttl = [];
end
opts_ttl = opts.ttl;

% assign defaults to the titles
if isempty(opts_ttl)
    opts_ttl = [1, 1.1, 1.2, 1.3];
end

titex = {strcat('(a) $\sigma_{{g,pp,ens|}_i}$ =', {' '}, num2str(opts_ttl(1),'%.1f')),...
    strcat('(b) $\sigma_{{g,pp,ens|}_i}$ =', {' '}, num2str(opts_ttl(2),'%.1f')),...
    strcat('(c) $\sigma_{{g,pp,ens|}_i}$ =', {' '}, num2str(opts_ttl(3),'%.1f')),...
    strcat('(d) $\sigma_{{g,pp,ens|}_i}$ =', {' '}, num2str(opts_ttl(4),'%.1f'))};

tbl = []; % initialize the datasheet to be plotted

n_dat = zeros(6,4); % number of data points in each boxplot

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

sigmas = zeros(4,1); % the ensemble average standard deviations obtained...
    % ...after post-flame lognormal sampling

% extract data from inputs
for j = 1 : 4
    dpp_ens = cell2mat(parsdata{j}(1).pp); 
    dpp_ens = dpp_ens(:,2);
    sigmas(j) = UTILS.GEOSTD(dpp_ens);
    
    for i = 1 : 6
        n_dat(i,j) = size(parsdata{j}(i).dpp_g, 1);
        tbl_temp = [repmat(t_id(7 - i), n_dat(i,j), 1),...
            repmat(sigmas(j), n_dat(i,j), 1), parsdata{j}(i).dpp_g(:,2)];
        tbl = vertcat(tbl, tbl_temp);
        
        if j == 1
            if i == 1
                legtxt{i} = strcat(' $n_{agg}/n_{agg_0}$ =', {' '},...
                    num2str(t_id(i), '%.0f'));
            else
                legtxt{i} = strcat(' $n_{agg}/n_{agg_0}$ =', {' '},...
                    num2str(t_id(i), '%.2f'));
            end
        end
    end
end

x_bc = categorical(round(tbl(:,2),2)); % generate the plot x axis

% make the boxplot
bc = boxchart(x_bc, tbl(:,3), 'GroupByColor', tbl(:,1));
for i = 1 : 6
    bc(i).BoxFaceColor = mc(i,:);
    bc(i).MarkerColor = mc(i,:);
    bc(i).MarkerSize = 5;
end
hold on

dx = 1/6; % gap to consider between the boxes to locate the scatter data

% put scatter data on top of boxplots
if ismember(opts_scat, {'ON', 'On', 'on'})
    for j = 1 : 4
        for i = 1 : 6
            x_scat = repmat(j + ((2 * i - 1) / 2 - 3) * dx, n_dat(i,j), 1);
            scatter(x_scat, parsdata{j}(i).dpp_g(:,2), 5, mc(i,:),...
                'MarkerEdgeAlpha', 0.3)
        end
    end
end

plot([1.5,1.5], [1,1.7], [2.5,2.5], [1,1.7], [3.5,3.5], [1,1.7],...
    'Color', 'k', 'LineWidth', 0.75)

box on
set(gca, 'TickLabelInterpreter','latex', 'FontSize', 18, 'TickLength', [0.005 0.005])
xlabel('$\sigma_{g,pp,ens}$ [-]', 'interpreter', 'latex', 'FontSize', 20)
ylabel('$\sigma_{g,pp,agg}$ [-]', 'interpreter', 'latex', 'FontSize', 20)
ylim([1, (ceil(max(ylim) / 0.1)) * 0.1 - 0.05])
text([0.6,1.6,2.6,3.6], repmat(1.61,4,1), titex, 'interpreter', 'latex', 'FontSize', 18)
legend(bc, cat(2, legtxt{:}), 'interpreter', 'latex', 'FontSize', 18,...
    'Orientation', 'horizontal', 'Location', 'northoutside', 'NumColumns', 3)

end

