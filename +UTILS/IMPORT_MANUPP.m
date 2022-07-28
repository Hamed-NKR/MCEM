function pp = IMPORT_MANUPP(fname, da, opts)
% "KINETIC" calculates the properties related to kinetics of aggrgeation.
% ----------------------------------------------------------------------- %
%
% Inputs:
%   fname: A cell array of name of the files to be imported
%   da: Projected area diameter of aggs
%   opts: Data import options
% ----------------------------------------------------------------------- %
%
% Output:
%   pp: Primary particle data structure
% ----------------------------------------------------------------------- %

n_agg = length(fname); % Number of aggs to be analyzed

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
    opts.visual = [];
end
opts_visual = opts.visual;
if isempty(opts_visual)
    opts_visual = 'off'; % Default not to plot the outputs
end

% Import the data and calculate properties
for i = 1 : n_agg
    f_ad = strcat('inputs\', fname{i}, '.csv'); % Making file address
    ppdat = readmatrix(f_ad); % Scanning data
    
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
    p1 = plot(da_uc, dpp_uc, 'Color', [0 0.4470 0.7410],...
        'LineStyle', '-.', 'LineWidth', 2.5);
    hold on
    
    % Plot manually sized aggs
    p2 = scatter(da, pp.mu_d(:,2), 25, [0.8500 0.3250 0.0980], '^');
    
    if (~isfield(opts, 'eb')) || isempty(opts.eb)
        opts.eb = 'on'; % default to plot errorbars
    end
    opts_eb = opts.eb;
    if ismember(opts_eb, {'ON', 'On', 'on'})
        e_p = pp.mu_d(:,2) .* abs(pp.std_d(:,2) - 1);
        e_n = pp.mu_d(:,2) .* abs(1 - 1 ./ pp.std_d(:,2));
        eb = errorbar(da, pp.mu_d(:,2), e_n, e_p, 'Marker', 'none',...
            'LineStyle', 'none');
        eb.Color = [0.8500 0.3250 0.0980];
    end
    
    if (~isfield(opts, 'label')) || isempty(opts.label)
        opts.label = 'on'; % default to plot std labels
    end
    opts_label = opts.label;
    if ismember(opts_label, {'ON', 'On', 'on'})
        text(da, pp.mu_d(:,2), num2str(pp.std_d(:,2), '%.2f'),...
            'VerticalAlignment', 'top', 'HorizontalAlignment', 'left',...
            'FontSize', 8)
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
    legend([p1, p2], {'Universal correlation', 'Manual'},...
        'Location', 'northwest', 'FontName', 'SansSerif', 'FontSize', 12);
    title('Primary particle size vs projected area equivalent size',...
        'FontName', 'SansSerif', 'FontWeight', 'bold', 'FontSize', 16)
end

end

