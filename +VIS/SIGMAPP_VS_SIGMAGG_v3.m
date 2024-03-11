function h = SIGMAPP_VS_SIGMAGG_v3(parsdata, opts)
% "SIGMAPP_VS_SIGMAGG_v3" plots boxcharts from temproal evolution of...
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
h.Position = [0, 0, 800, 600];
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
    opts_ttl = [1, 1.3];
end

titex = {strcat('(a) $\sigma_{{g,pp,ens|}_i}$ =', {' '}, num2str(opts_ttl(1),'%.1f')),...
    strcat('(b) $\sigma_{{g,pp,ens|}_i}$ =', {' '}, num2str(opts_ttl(2),'%.1f'))};

tbl = []; % initialize the datasheet to be plotted

n_dat = zeros(6,2); % number of data points in each boxplot

r_t = (0.02 / 0.5)^(1 / 4);
t_id = 0.5 * ones(5,1);
for i = 2 : 5
    t_id(i) = t_id(i) * r_t^(i-1);
end
t_id = [1; t_id]; % time identifier of data to be plotted 

% set box and point colors
mc = colormap(hot); % choose colormap
ii = round(1 + (length(mc) - 1) .* (0.05 : 0.7 / 5 : 0.75)'); % define the sample points over colormap range
mc = mc(ii,:); % sample across the colormap for different lifestages
% mc = flip(mc,1); % reverse the colormap direction (if needed)

sigmas = zeros(2,1); % the ensemble average standard deviations obtained...
    % ...after post-flame lognormal sampling

% extract data from inputs
for j = 1 : 2
    dpp_ens = cell2mat(parsdata{j}(1).pp); 
    dpp_ens = dpp_ens(:,2);
    sigmas(j) = UTILS.GEOSTD(dpp_ens);
    
    for i = 1 : 6
        n_dat(i,j) = size(parsdata{j}(i).dpp_g, 1);
        tbl_temp = [repmat(t_id(7 - i), n_dat(i,j), 1),...
            repmat(sigmas(j), n_dat(i,j), 1), parsdata{j}(i).dpp_g(:,2)];
        tbl = vertcat(tbl, tbl_temp);        
    end
end

x_bc = categorical(round(tbl(:,2),2)); % generate the plot x axis

% make the boxplot
bc = boxchart(x_bc, tbl(:,3), 'GroupByColor', tbl(:,1));
for i = 1 : 6
    if i < 6
        bc(i).BoxFaceColor = mc(i,:);
        bc(i).MarkerColor = mc(i,:);
    else
        bc(i).BoxFaceColor = [236,230,61] / 255;
        bc(i).MarkerColor = [236,230,61] / 255;        
    end
    bc(i).MarkerSize = 2.5;
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

plot([1.5,1.5], [0.95,1.7], 'Color', 'k', 'LineWidth', 0.5) % panel divider

% generate asymptotes
y_lim = unique(round(tbl(:,2),2));
p_lim = plot([0,1.5], [y_lim(1),y_lim(1)], 'Color', 'k', 'LineStyle', ':',...
    'LineWidth', 1.5);
plot([1.5,2.5], [y_lim(2),y_lim(2)], 'Color', 'k', 'LineStyle', ':',...
    'LineWidth', 1.5)

box on
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 18, 'TickLength', [0.01 0.01],...
    'YMinorTick', 'on', 'XTick', [])
xlabel({string(newline);string(newline);'$n_{agg}/n_{agg_0}$ [-]'}, 'interpreter', 'latex', 'FontSize', 20)
ylabel('$\overline{\sigma}_{g,pp,agg}$ [-]', 'interpreter', 'latex', 'FontSize', 20)
ylim([0.98, 1.62])
text([0.68,1.7], repmat(1.585,2,1), titex, 'interpreter', 'latex', 'FontSize', 20)

% xticks
plot([0.5833, 0.5833], [0.98, 0.9885], [0.75, 0.75], [0.98, 0.9885],[0.9167, 0.9167], [0.98, 0.9885],...
    [1.0833, 1.0833], [0.98, 0.9885], [1.25, 1.25], [0.98, 0.9885], [1.4167, 1.4167], [0.98, 0.9885],...
    'Color', 'k', 'LineWidth', 0.5)
plot([1.5833, 1.5833], [0.98, 0.9885], [1.75, 1.75], [0.98, 0.9885],[1.9167, 1.9167], [0.98, 0.9885],...
    [2.0833, 2.0833], [0.98, 0.9885], [2.25, 2.25], [0.98, 0.9885], [2.4167, 2.4167], [0.98, 0.9885],...
    'Color', 'k', 'LineWidth', 0.5)
plot([0.5833, 0.5833], [1.6115, 1.62], [0.75, 0.75], [1.6115, 1.62],[0.9167, 0.9167], [1.6115, 1.62],...
    [1.0833, 1.0833], [1.6115, 1.62], [1.25, 1.25], [1.6115, 1.62], [1.4167, 1.4167], [1.6115, 1.62],...
    'Color', 'k', 'LineWidth', 0.5)
plot([1.5833, 1.5833], [1.6115, 1.62], [1.75, 1.75], [1.6115, 1.62],[1.9167, 1.9167], [1.6115, 1.62],...
    [2.0833, 2.0833], [1.6115, 1.62], [2.25, 2.25], [1.6115, 1.62], [2.4167, 2.4167], [1.6115, 1.62],...
    'Color', 'k', 'LineWidth', 0.5)
text([0.565, 0.7, 0.85, 1.035, 1.18, 1.35, 1.565, 1.7, 1.85, 2.035, 2.18, 2.35],...
    repmat(0.957,12,1), {'1','0.5','0.22','0.1','0.04','0.02','1','0.5','0.22','0.1','0.04','0.02'},...
    'interpreter', 'latex', 'FontSize', 18)

legend(p_lim, '$\sigma_{g,pp,ens}$', 'interpreter', 'latex', 'FontSize', 20,...
    'Location', 'southeast')

end

