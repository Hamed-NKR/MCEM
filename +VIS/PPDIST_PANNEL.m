function h = PPDIST_PANNEL(pp_sim, ppdat0, Aggs, inds, fd0, opts)
% "PPDIST_PANNEL_SIM" plots a pannel of primary particle size...
%   ...distributions for simulated and experimental hybrid and non-hybrid...
%   ...aggregates along with their rendered strutcures.
% ----------------------------------------------------------------------- %
% 
% Inputs:
%   pp_sim: a cell array of primary particle datasheets for a uniform and...
%       ...a hybrid simulated aggregate.
%   ppdat0: File names of primary particle datasheets
%   Aggs: Structure containing aggregate data from TEM images (by ATEMS)
%   inds: Data indices to be analyzed in Aggs_im (in order of ppdat0)
%   fd0: Folder address for the primary particle data
% ----------------------------------------------------------------------- %
%
% Outputs:
%   h: output figure handle
% ----------------------------------------------------------------------- %

% initialize figure 
figure;
h = gcf;
h.Position = [0, 0, 1400, 750];
set(h, 'color', 'white');

% initialize layout
tt = tiledlayout(2,4);
tt.TileSpacing = 'tight';
tt.Padding = 'compact';

dist = cell(4,1); % placeholder for the aggregates' pp size distribution data

p = cell(2,4); % placeholder for subplots

% titles for subplots
subtitex = {'Uniform', 'Hybrid'};

% Set plot defaults
if ~exist('opts', 'var') 
    opts = struct();
end

% define number of bins for the size distribution histogram
if ~isfield(opts, 'n_bin')
    opts.n_bin = [];
end
n_bin = opts.n_bin;
if isempty(n_bin)
    n_bin = 10;
end

% Initialize data fields in the pp structure
pp_exp.mu_d = zeros(2,2); % Mean primary particle diameter of aggs (regular + geometric)
pp_exp.std_d = zeros(2,2); % Standard deviation of pp size (regular + geometric)
pp_exp.n = zeros(2,1); % Number of pps
pp_exp.d = cell(2,1); % Projected area diameter of pps
pp_exp.r = cell(2,1); % pp centroid locations
pp_exp.c = cell(2,1); % pp circularities

% Read pp data files
if ~exist('fd0', 'var') || isempty(fd0)
    fd0 = strcat(pwd, '\inputs');
end


for i = 1 : 4
    nexttile(i) % graphical structure demonstrations
    
    if i <= 2
        UTILS.PLOTPP(pp_sim{i}(:,3), pp_sim{i}(:,4), pp_sim{i}(:,5), pp_sim{i}(:,2),...
            [255, 210, 10] / 255)
        title(strcat(string(newline), subtitex{i}), 'FontSize', 20, 'interpreter','latex')
        
    else
        imshow(Aggs(inds(i-2)).image)
        axis('on', 'image')
        set(gca, 'xtick',[], 'xticklabel', [], 'ytick',[], 'yticklabel',[])

        title(strcat(string(newline), subtitex{i - 2}), 'FontSize', 20, 'interpreter','latex')        
    end
    
    nexttile(i + 4) % size distribution plots
    
    if i <= 2
        dist{i} = PAR.DISTPP(pp_sim{i}(:,2), n_bin); % get the distibution from pp info matrix
        
    else
        fd = strcat(fd0, '\', ppdat0{i-2}, '.csv'); % make the file address
        ppdat = readmatrix(fd); % scan the data
        
        % remove the header and the scale detector and save other data
        pp_exp.n(i) = size(ppdat, 1) - 2;
        pp_exp.d{i} = 2 * sqrt(ppdat(3:end, 2) / pi);
        pp_exp.r{i} = ppdat(3:end, 3:4);
        pp_exp.c{i} = 2 * sqrt(pi * ppdat(3:end, 2)) ./ ppdat(3:end, 5);

        % calculate mean and std of size
        pp_exp.mu_d(i,:) = [mean(pp_exp.d{i}), geomean(pp_exp.d{i})];
        pp_exp.std_d(i,:) = [std(pp_exp.d{i}), UTILS.GEOSTD(pp_exp.d{i})];

        pp_exp.fname{i} = ppdat0{i - 2}; % save file name for future reference
        
        dist{i} = PAR.DISTPP(pp_exp.d{i}, n_bin); % get the distibution from pp info matrix        
    end

    % plot the distribution histogram
    xhist = [dist{i}.d_bin(1); repelem(dist{i}.d_bin(2 : end - 1), 2);...
        dist{i}.d_bin(end)];
    if i <= 2; xhist = 1e9 * xhist; end
    yhist = repelem(dist{i}.dn_dlogdpp, 2);
    p{i,2} = plot(xhist, yhist, 'Color', [0.1 0.1 0.1], 'LineWidth', 1);
    
    % set subplot axes
    box on
    set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 18,...
        'TickLength', [0.02 0.02], 'XScale', 'log')
    xlim([min(xhist), max(xhist)])
    ylim([0, ceil(max(yhist) / 100) * 100 + 50])
    if i == 1
        ylabel('d$n_{pp}$/dlog($d_{pp}$) [-]', 'interpreter', 'latex',...
            'FontSize', 20)
    else
        set(gca, 'yticklabel',[])
    end
    
end

xlabel(tt, '$d_{pp}$ [nm]', 'interpreter', 'latex', 'FontSize', 20)
title(tt, '$\ \ $Simulated$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ $Experimental',...
    'FontSize', 24, 'interpreter','latex')

end

