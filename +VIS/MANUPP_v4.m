function pp = MANUPP_v4(fname_pp, da, fad_pp, n_hyb, opts)
% "manu_pp" calculates the properties of primary particles extracted from...
%   ...manual sizing of a real TEM aggregate.
% ----------------------------------------------------------------------- %
%
% Inputs:
%   fname_pp: A cell array of name of the files to be imported
%   da: Projected area diameter of aggs
%   fad_pp: folder address (e.g., 'inputs\manual\')
%   n_hyb: Number of hybridity regions
%   opts: Data import options
% ----------------------------------------------------------------------- %
%
% Output:
%   pp: Primary particle data structure
% ----------------------------------------------------------------------- %

n_agg = length(fname_pp); % Number of aggs to be analyzed

% Initialize the pp structure
pp = struct();
pp.mu_d = zeros(n_agg,2); % Mean primary particle diameter of aggs (regular + geometric)
pp.std_d = zeros(n_agg,2); % Standard deviation of pp size (regular + geometric)
pp.n = zeros(n_agg,1); % Number of pps
pp.d = cell(n_agg,1); % Projected area diameter of pps
pp.r = cell(n_agg,1); % pp centroid locations
pp.c = cell(n_agg,1); % pp circularities

% Set defaults for dp vs da plotting
if ~exist('opts', 'var') 
    opts = struct();
end

if ~isfield(opts, 'visual')
    opts.visual = []; % define plotting decision variable
end
opts_visual = opts.visual;
if isempty(opts_visual)
    opts_visual = 'on'; % default to plot the outputs
end

if (~exist('n_hyb', 'var')) || isempty(n_hyb)
    n_hyb = ones(n_agg,1); % default for all to be non-hybrid
end

p = cell(4,1); % placeholder for plots
legtxt = cell(4,1); % legend text placeholder

pp.fname = fname_pp; % store file names

ii = cell(2,2); % second-hand index for source-based filtering of scattered data

% Import the data and calculate properties
for i = 1 : n_agg
    fad = strcat(fad_pp, '\', fname_pp{i}, '.csv'); % Making file address
    ppdat = readmatrix(fad); % Scanning data
    
    pp.n(i) = size(ppdat, 1) - 2; % Removing the header and the scale detector 
    pp.d{i} = 2 * sqrt(ppdat(3:end, 2) / pi);
    pp.r{i} = ppdat(3:end, 3:4);
    pp.c{i} = 2 * sqrt(pi * ppdat(3:end, 2)) ./ ppdat(3:end, 5);
    
    pp.mu_d(i,:) = [mean(pp.d{i}), geomean(pp.d{i})];
    pp.std_d(i,:) = [std(pp.d{i}), UTILS.GEOSTD(pp.d{i})];
    
    % detect source and hybridity and label aggregates
    if (contains(fname_pp{i},'D1') || contains(fname_pp{i},'D2') ||...
        contains(fname_pp{i},'D3')) && (n_hyb(i) == 1) % Argonaut, semiuniform
        ii{1,1} = [ii{1,1}, i];
    
    elseif (contains(fname_pp{i},'D1') || contains(fname_pp{i},'D2') ||...
        contains(fname_pp{i},'D3')) && (n_hyb(i) > 1) % Argonaut, biregional hybrid
        ii{1,2} = [ii{1,2}, i];
    
%     elseif (contains(fname_pp{i},'D1') || contains(fname_pp{i},'D2') ||...
%         contains(fname_pp{i},'D3')) && (n_hyb(i) > 2) % Argonaut, multiregional hybrid
%         ii{1,3} = [ii{1,3}, i];
%     
    elseif contains(fname_pp{i},'D8') && (n_hyb(i) == 1) % Flare, semiuniform
        ii{2,1} = [ii{2,1}, i];
    
    elseif contains(fname_pp{i},'D8') && (n_hyb(i) > 1) % Flare, biregional hybrid
        ii{2,2} = [ii{2,2}, i];

%     elseif (contains(fname_pp{i},'D8')) && (n_hyb(i) > 2) % Flare, multiregional hybrid
%         ii{2,3} = [ii{2,3}, i];    
    end
end

if ismember(opts_visual, {'ON', 'On', 'on'})
    % initialize figure properties
    figure
    h = gcf;
    if ~all(h.Position == [0, 0, 850, 500])
        h.Position = [0, 0, 850, 500];
    end
    set(h, 'color', 'white')
    
    % initialize layout
    tt = tiledlayout(1,4);
    tt.TileSpacing = 'compact';
    tt.Padding = 'compact';
    nexttile(2,[1,2])
    
    % Plot universal correlation
    dpp_uc = linspace(5, 65, 1000);
    da_uc = 100 * (dpp_uc / 17.8).^(1 / 0.35);
    p0 = plot(da_uc, dpp_uc, 'Color', [0.4940 0.1840 0.5560],...
        'LineStyle', '-.', 'LineWidth', 3);
    hold on
    
    % plot manual dpp vs da data
    p{1} = scatter(da(ii{1,1}), pp.mu_d(ii{1,1},2), 60, [0.1 0.1 0.1], 's');
    legtxt{1} = 'Argonaut - Semi-uniform';
    
    p{2} = scatter(da(ii{1,2}), pp.mu_d(ii{1,2},2), 80, [0.1 0.1 0.1], 's', 'filled');
    legtxt{2} = 'Argonaut - Hybrid';
    
    p{3} = scatter(da(ii{2,1}), pp.mu_d(ii{2,1},2), 40, [0.1 0.1 0.1], '^');
    legtxt{3} = 'Flare - Semi-uniform';
    
    p{4} = scatter(da(ii{2,2}), pp.mu_d(ii{2,2},2), 50, [0.1 0.1 0.1], '^', 'filled');
    legtxt{4} = 'Flare - Hybrid';
    
    if (~isfield(opts, 'select')) || isempty(opts.select)
        opts.select = 'off'; % default not to select points
    end
    opts_select = opts.select;
    if ismember(opts_select, {'ON', 'On', 'on'})
        if (~isfield(opts, 'tem')) || isempty(opts.tem)
            opts.tem = {'D1_062_1_S', 'D3_044', 'D3_072', 'D8_026',...
                'D8_052_1_S_H', 'D8_051', 'D8_052_2_L_U'};
        end
        opts_tem = opts.tem;
        
        jj = zeros(length(opts_tem),1);
        for j = 1 : length(opts_tem)
            for i = 1 : n_agg
                if contains(fname_pp{i}, opts_tem{j})
                    jj(j) = i;
                end
            end
        end
        
        scatter(da(jj), pp.mu_d(jj,2), 150, [0.1 0.1 0.1], 'o',...
            'LineWidth', 1);
    end    
    
    if (~isfield(opts, 'label')) || isempty(opts.label)
        opts.label = 'off'; % default not to plot std labels
    end
    opts_label = opts.label;
    if ismember(opts_label, {'ON', 'On', 'on'})
        labtex = cell(n_agg,1);
        for i = 1 : n_agg
            labtex{i} = strrep(fname_pp{i}(4:end), '_', '-');
        end
        
        text(da, pp.mu_d(:,2), labtex, 'Interpreter', 'latex',...
            'VerticalAlignment', 'middle', 'HorizontalAlignment', 'center',...
            'FontSize', 8)        
    end        

    box on
    
    set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 16, 'TickLength', [0.02 0.02])
    
    xlabel('$d_a$ [nm]', 'Interpreter', 'latex', 'FontWeight', 'bold',...
        'FontSize', 20)
    xlim([30, 1300])
    set(gca, 'XScale', 'log')
    
    ylabel('$\overline{d}_{pp}$ [nm]', 'Interpreter', 'latex',...
        'FontWeight', 'bold', 'FontSize', 20)
    ylim([7, 70])
    set(gca, 'YScale', 'log')
    xticks([100, 1000])
    yticks(10 : 10 : 70)
    
    lgd = legend([p0, cat(1, p{:})'], cat(2, {'Universal correlation'}, legtxt{:}),...
        'interpreter', 'latex', 'Location', 'southoutside', 'FontSize', 16,...
        'Orientation', 'horizontal', 'NumColumns', 3);
end

% axis square
lgd.Layout.Tile = 'south';

end
