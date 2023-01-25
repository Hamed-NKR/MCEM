function f1 = DPP_VS_DA_v5(parsdata, opts)
% DPP_VS_DA_v5 plots of primary particle size vs. aggregate size as a...
%   ...function of time and number of clusters.
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

% initialize temporal figure 
f1 = figure(1);
f1.Position = [0, 0, 700, 700];
set(f1, 'color', 'white')

% initialize plot & legend placholders
p1 = cell(7,1);
legtxt1 = cell(7,1);
legtxt1{7} = 'Olfert $\&$ Rogak (2019)';

% set plotting options variable
if ~exist('opts', 'var') 
    opts = struct();
end

% reproduce saving criteria for legends
r_kk1 = (0.02 / 0.5)^(1 / 4);
kk1 = 0.5 * ones(5,1);
for i = 2 : 5
    kk1(i) = kk1(i) * r_kk1^(i-1);
end
kk1 = [1; kk1];

% make the colormap and initialize the maker properties
mc1 = colormap(hot);
ii1 = round(1 + (length(mc1) - 1) .* (0.05 : 0.7 / 5 : 0.75)');
mc1 = mc1(ii1,:);
ms = [20, 20, 20, 30, 25, 35];
mt = {'o', '^', 'v', 's', 'd', 'p'};
mc1(6,:) = [236,230,61] / 255;

% generate universal correlation data
dpp_uc = logspace(0, 4, 1000);
da_uc = 100 * (dpp_uc / 17.8).^(1 / 0.35);

p1{7,1} = plot(da_uc, dpp_uc, 'Color', [0.4940 0.1840 0.5560],...
    'LineStyle', '-.', 'LineWidth', 3);
hold on

for i = 1 : 6
    p1{i,1} = scatter(1e9 * parsdata(i).da, 1e9 * parsdata(i).dpp_g(:,1),...
        ms(i), mc1(i,:), mt{i}, 'LineWidth', 1); % plot temporal dpp vs. da
    
    % make legends for pannel 1
    if i == 1
        legtxt1(i) = strcat('$n_{agg}$/$n_{agg_0}$ =',...
            {' '}, num2str(kk1(i), '%.0f'));
    else
        legtxt1(i) = strcat('$n_{agg}$/$n_{agg_0}$ =',...
            {' '}, num2str(kk1(i), '%.2f'));
    end
end

box on
xlim([35, 2000])
ylim([10, 40])
ax1 = gca;
set(ax1, 'FontSize', 18, 'TickLength', [0.02 0.02], 'XScale', 'log',...
    'YScale', 'log', 'TickLabelInterpreter','latex')
ax1.XTick = [100, 1000];
ax1.YTick = [10,20, 30, 40, 50];

xlabel('$\overline{d}_a$ [nm]', 'FontSize', 20, 'interpreter','latex')
ylabel('$\overline{d}_{pp}$ [nm]', 'FontSize', 20, 'interpreter', 'latex')

legend(cat(1, p1{:})', cat(2,legtxt1(:)), 'Location', 'southeast',...
    'FontSize', 16, 'interpreter', 'latex');

% initialize figure 
f2 = figure(2);
f2.Position = [0, 0, 1000, 650];
set(f2, 'color', 'white')

% define plotting variables
p2 = cell(10,1);
legtxt2 = cell(10,1);
legtxt2{10} = 'Olfert $\&$ Rogak (2019)';
kk2 = [1, 2, 3, 4, 5, 7, 10, 15, inf];
mc2 = colormap(turbo);
ii2 = round(1 + (length(mc2) - 1) .* (0.05 : 0.9 / 7 : 0.95)');
mc2 = mc2(ii2,:);
mc2 = flip (mc2,1);
ms2 = [20, 20, 20, 30, 25, 35, 35, 30];
mt2 = {'o', '^', 'v', 's', 'd', 'p', 'h', '*'};

% concatinate the data
dpp_ens = cat(1,parsdata.dpp_g);
dpp_ens = dpp_ens(:,1);
nagg_ens = length(dpp_ens);
da_ens = cat(1,parsdata.da);
nhyb_ens = cat(1,parsdata.n_hyb);
pp_ens = cat(1,parsdata.pp);

% find and remove duplicates
ij = nchoosek(1 : nagg_ens, 2);
ind_flt = zeros(nagg_ens,1);
for i = 1 : length(ij)
    if isequal(sort(unique(pp_ens{ij(i,1)}(:,1))),...
            sort(unique(pp_ens{ij(i,2)}(:,1))))
        ind_flt(ij(i,2)) = 1;
    end
end
ind_flt = logical(ind_flt);
dpp_flt = dpp_ens(~ind_flt);
da_flt = da_ens(~ind_flt);
nhyb_flt = nhyb_ens(~ind_flt);
pp_flt = pp_ens(~ind_flt);

p2{10,1} = plot(da_uc, dpp_uc, 'Color', [0.4940 0.1840 0.5560],...
    'LineStyle', '-.', 'LineWidth', 3); % universal correlation
hold on

% initial ensemble average pp mean diameter
pp1_ens = cell2mat(parsdata(1).pp);
dpp1_ens = geomean(pp1_ens(:,2));
p2{9,1} = plot(da_uc, 1e9 * repelem(dpp1_ens, length(da_uc)),...
    'Color', [0.5, 0.5, 0.5], 'LineStyle', '--', 'LineWidth', 2);
legtxt2{9} = '$\overline{d}_{pp} = \overline{d}_{{pp,ens|}_0}$';
% legtxt2(9) = strcat('$\overline{d}_{{pp,ens|}_0}$ =', {' '},...
%     num2str(1e9 * dpp1_ens,'%.1f'), ' nm');

iii = cell(8,1); % index sorting placeholder based on number of internal clusters

for i = 1 : 8
    iii{i} = (nhyb_flt >= kk2(i)) & (nhyb_flt < kk2(i+1));
    
    p2{i,1} = scatter(1e9 * da_flt(iii{i}), 1e9 * dpp_flt(iii{i}),...
        ms2(i), mc2(i,:), mt2{i}, 'LineWidth', 1); % plot dpp vs. da as a function of internal cluster counts
    
    % make legends for pannel 2
    switch i
        case {1,2,3,4}
            legtxt2(i) = strcat('$n_{hyb}$ =',...
                {' '}, num2str(kk2(i), '%d'));
        case {5,6,7}
            legtxt2(i) = strcat(num2str(kk2(i), '%d'), {' '},...
                '$\leq n_{hyb} <$', {' '}, num2str(kk2(i+1), '%d'));
        otherwise
            legtxt2(i) = strcat('$n_{hyb} \geq$', {' '},...
                num2str(kk2(i), '%d'));
    end
end

box on
xlim([35, 2000])
ylim([10, 40])
ax2 = gca;
set(ax2, 'FontSize', 18, 'TickLength', [0.02 0.02], 'XScale', 'log',...
    'YScale', 'log', 'TickLabelInterpreter','latex')
ax2.XTick = [100, 1000];
ax2.YTick = [10, 20, 30, 40, 50];

xlabel('$\overline{d}_a$ [nm]', 'FontSize', 20, 'interpreter','latex')
ylabel('$\overline{d}_{pp}$ [nm]', 'FontSize', 20, 'interpreter', 'latex')

lgd2 = legend(cat(1, p2{1:10})', cat(2,legtxt2(1:10)), 'Location', 'northeastoutside',...
    'FontSize', 16, 'interpreter', 'latex');

% disable aggregate rendering option if not defined
if ~isfield(opts, 'render')
    opts.render = 'off';
end
opts_render = opts.render;

% identify and save aggregates of interest to be rendered
if ismember(opts_render, {'ON', 'On', 'on'})
    if ~isfield(opts, 'id') % set ids of aggs to be saved if not defined
        da1 = sort(da_flt(iii{1}));
        da3 = sort(da_flt(iii{3}));
        da5 = sort(da_flt(iii{5}));
        da8 = sort(da_flt(iii{8}));
        
        opts.id = [find(da_flt == da1(ceil(end/12)),1);...
            find(da_flt == da1(ceil(end/2)),1);...
            find(da_flt == da1(ceil(11*end/12)),1);...
            find(da_flt == da3(ceil(end/3)),1);...
            find(da_flt == da5(ceil(end/2)),1);...
            find(da_flt == da8(ceil(2*end/3)),1)];
    end
    opts_id = opts.id;
    
    % highlight the points of interest
    scatter(1e9 * da_flt(opts_id), 1e9 * dpp_flt(opts_id),...
        100, [0.1, 0.1, 0.1], 'o', 'LineWidth', 1);
%     text(x_annot(1:2), y_annot(1:2), '\downarrow', 'VerticalAlignment' , 'bottom',...
%         'FontSize', 16, 'FontWeight', 'bold')
%     text(x_annot(3:5), y_annot(3:5), '\uparrow', 'VerticalAlignment', 'cap',...
%         'FontSize', 16, 'FontWeight', 'bold')
    lgd2.String(11) = [];

    exportgraphics(f2, 'outputs\dpp-vs-da.emf', 'BackgroundColor','none',...
        'ContentType','vector')
    
    % render and save selected aggregates
    opts2.cc = 'on';
    opts2.cm = bone;
    
    f3 = figure(3);
    UTILS.PLOTPP(pp_flt{opts_id(1)}(:,3), pp_flt{opts_id(1)}(:,4),...
        pp_flt{opts_id(1)}(:,5), pp_flt{opts_id(1)}(:,2),...
        pp_flt{opts_id(1)}(:,6), opts2)
    exportgraphics(f3, 'outputs\render1-1.png', 'Resolution', 1200)

    f4 = figure(4);
    UTILS.PLOTPP(pp_flt{opts_id(2)}(:,3), pp_flt{opts_id(2)}(:,4),...
        pp_flt{opts_id(2)}(:,5), pp_flt{opts_id(2)}(:,2),...
        pp_flt{opts_id(2)}(:,6), opts2)
    exportgraphics(f4, 'outputs\render1-2.png', 'Resolution', 1200)

    f5 = figure(5);
    UTILS.PLOTPP(pp_flt{opts_id(3)}(:,3), pp_flt{opts_id(3)}(:,4),...
        pp_flt{opts_id(3)}(:,5), pp_flt{opts_id(3)}(:,2),...
        pp_flt{opts_id(3)}(:,6), opts2)
    exportgraphics(f5, 'outputs\render1-3.png', 'Resolution', 1200)

    f6 = figure(6);
    UTILS.PLOTPP(pp_flt{opts_id(4)}(:,3), pp_flt{opts_id(4)}(:,4),...
        pp_flt{opts_id(4)}(:,5), pp_flt{opts_id(4)}(:,2),...
        pp_flt{opts_id(4)}(:,6), opts2)
    exportgraphics(f6, 'outputs\render3.png', 'Resolution', 1200)

    f7 = figure(7);
    UTILS.PLOTPP(pp_flt{opts_id(5)}(:,3), pp_flt{opts_id(5)}(:,4),...
        pp_flt{opts_id(5)}(:,5), pp_flt{opts_id(5)}(:,2),...
        pp_flt{opts_id(5)}(:,6), opts2)
    exportgraphics(f7, 'outputs\render5.png', 'Resolution', 1200)

    f8 = figure(8);
    UTILS.PLOTPP(pp_flt{opts_id(6)}(:,3), pp_flt{opts_id(6)}(:,4),...
        pp_flt{opts_id(6)}(:,5), pp_flt{opts_id(6)}(:,2),...
        pp_flt{opts_id(6)}(:,6), opts2)
    exportgraphics(f8, 'outputs\render15.png', 'Resolution', 1200)
end

end

