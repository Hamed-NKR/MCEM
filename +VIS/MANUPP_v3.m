function pp = MANUPP_v3(fname_pp, da, fad_pp, n_hyb, opts)
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
    opts_visual = 'on'; % default not to plot the outputs
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
        contains(fname_pp{i},'D3')) && (n_hyb(i) == 1) % Argonaut, uniform
        ii{1,1} = [ii{1,1}, i];
    
    elseif (contains(fname_pp{i},'D1') || contains(fname_pp{i},'D2') ||...
        contains(fname_pp{i},'D3')) && (n_hyb(i) > 1) % Argonaut, hybrid
        ii{1,2} = [ii{1,2}, i];
    
    elseif contains(fname_pp{i},'D8') && (n_hyb(i) == 1) % Flare, uniform
        ii{2,1} = [ii{2,1}, i];
    
    elseif contains(fname_pp{i},'D8') && (n_hyb(i) > 1) % Flare, hybrid
        ii{2,2} = [ii{2,2}, i];
    end
end

if ismember(opts_visual, {'ON', 'On', 'on'})
    % initialize figure properties
    figure
    h = gcf;
    if ~all(h.Position == [0, 0, 700, 700])
        h.Position = [0, 0, 700, 700];
    end
    set(h, 'color', 'white')
    
    % Plot universal correlation
    dpp_uc = linspace(5, 65, 1000);
    da_uc = 100 * (dpp_uc / 17.8).^(1 / 0.35);
    p0 = plot(da_uc, dpp_uc, 'Color', [0.4940 0.1840 0.5560],...
        'LineStyle', '-.', 'LineWidth', 3);
    hold on
    
    p{1} = scatter(da(ii{1,1}), pp.mu_d(ii{1,1}), 60, pp.std_d(ii{1,1},2), 's');
    legtxt{1} = 'Argonaut burner - Uniform';
    
    p{2} = scatter(da(ii{1,2}), pp.mu_d(ii{1,2}), 60, pp.std_d(ii{1,2},2), 's', 'filled');
    legtxt{2} = 'Argonaut burner - Hybrid';
    
    p{3} = scatter(da(ii{2,1}), pp.mu_d(ii{2,1}), 40, pp.std_d(ii{2,1},2), '^');
    legtxt{3} = 'Flare - Uniform';
    
    p{4} = scatter(da(ii{2,2}), pp.mu_d(ii{2,2}), 40, pp.std_d(ii{2,2},2), '^', 'filled');
    legtxt{4} = 'Flare - Hybrid';
    
    cb = colorbar;
    colormap turbo
    cb.Label.String = '$\sigma_{g,pp,agg}$ [-]';
    cb.Label.Interpreter  = 'latex';
    cb.Label.Rotation = 0;
    cb.TickLabelInterpreter  = 'latex';
    cb.FontSize = 14;
    cb.Label.FontSize = 18;
    cbpos = get(cb, 'Position');
    cb.Label.Position = [cbpos(1) + 1.3 , cbpos(2) + 1.61];
        
    if (~isfield(opts, 'label')) || isempty(opts.label)
        opts.label = 'off'; % default to plot std labels
    end
    opts_label = opts.label;
    if ismember(opts_label, {'ON', 'On', 'on'})
        labtex = cell(n_agg,1);
        for i = 1 : n_agg
            labtex{i} = strcat(fname_pp{i}(4:5), '-', fname_pp{i}(7:9));
        end
        
        text(da, pp.mu_d(:,2), labtex, 'VerticalAlignment', 'bottom',...
            'HorizontalAlignment', 'right', 'FontSize', 8)        
    end        

    box on
    
    set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 16, 'TickLength', [0.01 0.01])
    
    xlabel('$d_a$ [nm]', 'Interpreter', 'latex', 'FontWeight', 'bold',...
        'FontSize', 20)
    xlim([30, 1300])
    set(gca, 'XScale', 'log')
    
    ylabel('$\overline{d}_{pp}$ [nm]', 'Interpreter', 'latex',...
        'FontWeight', 'bold', 'FontSize', 20)
    ylim([7, 75])
    set(gca, 'YScale', 'log')
    yticks(10 : 10 : 70)
    
    legend([cat(1, p{:})', p0], cat(2, legtxt{:}, {'Universal correlation'}),...
        'interpreter', 'latex', 'Location', 'northwest', 'FontSize', 14)
%         'Orientation', 'horizontal', 'NumColumns', 2)
end

end
