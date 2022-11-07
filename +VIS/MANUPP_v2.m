function pp = MANUPP_v2(fname0, da, fad0, n_hyb, opts)
% "manu_pp" calculates the properties of primary particles extracted from...
%   ...manual sizing of a real TEM aggregate.
% ----------------------------------------------------------------------- %
%
% Inputs:
%   fnam0: A cell array of name of the files to be imported
%   da: Projected area diameter of aggs
%   fad0: folder address (e.g., 'inputs\manual\')
%   n_hyb: Number of hybridity regions
%   opts: Data import options
% ----------------------------------------------------------------------- %
%
% Output:
%   pp: Primary particle data structure
% ----------------------------------------------------------------------- %

n_agg = length(fname0); % Number of aggs to be analyzed

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
    opts_visual = 'on'; % default not to plot the outputs
end

if ~isfield(opts, 'lbl')
    opts.lbl = cell(1, n_agg); % define aggregate labeling parameter
end

if (~exist('n_hyb', 'var')) || isempty(n_hyb)
    n_hyb = ones(n_agg,1); % default for all to be non-hybrid
end

p = cell(3,1); % placeholder for plots
legtxt = cell(3,1); % legend text placeholder

pp.fname = fname0; % store file names

% Import the data and calculate properties
for i = 1 : n_agg
    fad = strcat(fad0, '\', fname0{i}, '.csv'); % Making file address
    ppdat = readmatrix(fad); % Scanning data
    
    pp.n(i) = size(ppdat, 1) - 2; % Removing the header and the scale detector 
    pp.d{i} = 2 * sqrt(ppdat(3:end, 2) / pi);
    pp.r{i} = ppdat(3:end, 3:4);
    pp.c{i} = 2 * sqrt(pi * ppdat(3:end, 2)) ./ ppdat(3:end, 5);
    
    pp.mu_d(i,:) = [mean(pp.d{i}), geomean(pp.d{i})];
    pp.std_d(i,:) = [std(pp.d{i}), UTILS.GEOSTD(pp.d{i})];
end

if ismember(opts_visual, {'ON', 'On', 'on'})
    % initialize figure properties
    figure
    h = gcf;
    if ~all(h.Position == [0, 0, 600, 600])
        h.Position = [0, 0, 600, 600];
    end
    set(h, 'color', 'white')
    
    % Plot universal correlation
    dpp_uc = linspace(5, 65, 1000);
    da_uc = 100 * (dpp_uc / 17.8).^(1 / 0.35);
    p1 = plot(da_uc, dpp_uc, 'Color', [0.5 0.5 0.5],...
        'LineStyle', '-.', 'LineWidth', 2.5);
    hold on
    
    p{1} = scatter(da(n_hyb == 1), pp.mu_d(n_hyb == 1), 35, [0.8500 0.3250 0.0980], '^');
    legtxt{1} = 'Near-uniform';
    
    p{2} = scatter(da(n_hyb == 2), pp.mu_d(n_hyb == 2), 50, [0.3010 0.7450 0.9330], 's');
    legtxt{2} = 'Hybrid, \itn_{hyb,vis}\rm = 2';
    
    p{3} = scatter(da(n_hyb > 2), pp.mu_d(n_hyb > 2), 50, [0 0.4470 0.7410], 'o');
    legtxt{3} = 'Hybrid, \itn_{hyb,vis}\rm > 2';

    if (~isfield(opts, 'eb')) || isempty(opts.eb)
        opts.eb = 'on'; % default to plot errorbars
    end
    opts_eb = opts.eb;
    if ismember(opts_eb, {'ON', 'On', 'on'})
        % calculate error bounds
        e_p = pp.mu_d(:,2) .* abs(pp.std_d(:,2) - 1);
        e_n = pp.mu_d(:,2) .* abs(1 - 1 ./ pp.std_d(:,2));
        
        eb = cell(3,1);
        
        eb{1} = errorbar(da(n_hyb == 1), pp.mu_d(n_hyb == 1),...
            e_n(n_hyb == 1), e_p(n_hyb == 1), 'Marker', 'none',...
            'LineStyle', 'none');
        eb{1}.Color = [0.8500 0.3250 0.0980];
        
        eb{2} = errorbar(da(n_hyb == 2), pp.mu_d(n_hyb == 2),...
            e_n(n_hyb == 2), e_p(n_hyb == 2), 'Marker', 'none',...
            'LineStyle', 'none');
        eb{2}.Color = [0.3010 0.7450 0.9330];

        eb{3} = errorbar(da(n_hyb > 2), pp.mu_d(n_hyb > 2),...
            e_n(n_hyb > 2), e_p(n_hyb > 2), 'Marker', 'none',...
            'LineStyle', 'none');
        eb{3}.Color = [0 0.4470 0.7410];        
    end
    
    if (~isfield(opts, 'label')) || isempty(opts.label)
        opts.label = 'on'; % default to plot std labels
    end
    opts_label = opts.label;
    if ismember(opts_label, {'ON', 'On', 'on'})
        text(da, pp.mu_d(:,2), num2str(pp.std_d(:,2), '%.2f'),...
            'VerticalAlignment', 'top', 'HorizontalAlignment', 'left',...
            'FontSize', 8)
        text(da, pp.mu_d(:,2), opts.lbl, 'VerticalAlignment', 'bottom',...
            'HorizontalAlignment', 'right', 'FontSize', 8)        
    end        

    box on
    set(gca, 'FontName', 'SansSerif', 'FontSize', 12, 'TickLength', [0.02 0.02])
    xlabel('d_a (nm)', 'FontName', 'SansSerif', 'FontWeight', 'bold',...
        'FontSize', 14)
    ylabel('$\overline {d_{pp}} (nm)$', 'Interpreter', 'latex',...
        'FontName', 'SansSerif', 'FontWeight', 'bold', 'FontSize', 14)
    xlim([5, inf])
    ylim([5, 65])
    set(gca, 'XScale', 'log')
    set(gca, 'YScale', 'log')
    legend([p1, cat(1, p{:})'], cat(2, {'Universal correlation'}, legtxt{:}),...
        'Location', 'northwest', 'FontName', 'SansSerif', 'FontSize', 12);
    title('Primary particle size vs projected area equivalent size',...
        'FontName', 'SansSerif', 'FontWeight', 'bold', 'FontSize', 16)
end

end
