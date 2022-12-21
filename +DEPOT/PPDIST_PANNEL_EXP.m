function h = PPDIST_PANNEL_EXP(ppdat0, Aggs, inds, fd0)
% "PPDIST_PANNEL_EXP" plots a pannel of primary particle size...
%   ...distributions for experimental hybrid and non-hybrid aggregates...
%   ...along with their TEM images.
% ----------------------------------------------------------------------- %
% 
% Inputs:
%   ppdat0: File names of primary particle datasheets
%   Aggs: Structure containing aggregate data from TEM images (by ATEMS)
%   inds: Data indices to be analyzed in Aggs_im (in order of ppdat0)
%   fd0: Folder address for the primary particle data
% ----------------------------------------------------------------------- %
%
% Outputs:
%   h: Output figure handle
% ----------------------------------------------------------------------- %

% initialize figure 
figure;
h = gcf;
h.Position = [0, 0, 800, 700];
set(h, 'color', 'white');

% initialize layout
tt = tiledlayout(2,2);
tt.TileSpacing = 'compact';
tt.Padding = 'compact';

% Initialize the pp structure
pp = struct();
pp.mu_d = zeros(2,2); % Mean primary particle diameter of aggs (regular + geometric)
pp.std_d = zeros(2,2); % Standard deviation of pp size (regular + geometric)
pp.n = zeros(2,1); % Number of pps
pp.d = cell(2,1); % Projected area diameter of pps
pp.r = cell(2,1); % pp centroid locations
pp.c = cell(2,1); % pp circularities

% Read pp data files
if ~exist('fd0', 'var') || isempty(fd0)
    fd0 = strcat(pwd, '\inputs');
end
for i = 1 : 2
    fd = strcat(fd0, '\', ppdat0{i}, '.csv'); % make the file address
    ppdat = readmatrix(fd); % scan the data
    
    % remove the header and the scale detector and save other data
    pp.n(i) = size(ppdat, 1) - 2;
    pp.d{i} = 2 * sqrt(ppdat(3:end, 2) / pi);
    pp.r{i} = ppdat(3:end, 3:4);
    pp.c{i} = 2 * sqrt(pi * ppdat(3:end, 2)) ./ ppdat(3:end, 5);
    
    % calculate mean and std of size
    pp.mu_d(i,:) = [mean(pp.d{i}), geomean(pp.d{i})];
    pp.std_d(i,:) = [std(pp.d{i}), UTILS.GEOSTD(pp.d{i})];
    
    pp.fname{i} = ppdat0{i}; % save file name for future reference
end

dist = cell(2,1); % placeholder for the aggregates' pp size distribution data

p = cell(2, 2); % placeholder for subplots

% titles for subplots
titex = {'Experimental - Uniform', 'Experimental - Hybrid'};

for i = 1 : 2
    nexttile(2 + i) % distribution plot
    
    dist{i} = DISTPP(pp.d{i}, 10); % get the distibution from pp info matrix
    
    % plot the distribution histogram
    xhist = [dist{i}.d_bin(1); repelem(dist{i}.d_bin(2 : end - 1), 2);...
        dist{i}.d_bin(end)];
    yhist = repelem(dist{i}.dn_dlogdpp, 2);
    p{i,2} = plot(xhist, yhist, 'Color', [0.1 0.1 0.1], 'LineWidth', 1);
    
    % set subplot axes
    box on
    set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 18,...
        'TickLength', [0.02 0.02], 'XScale', 'log')
    xlim([min(xhist), max(xhist)])
    xlabel('$d_{pp}$ [nm]', 'interpreter', 'latex', 'FontSize', 20)
    ylim([0,950])
    if i == 1
        ylabel('d$n_{pp}$/dlog($d_{pp}$) [-]', 'interpreter', 'latex',...
            'FontSize', 20)
    else
        set(gca, 'yticklabel',[])
    end
    
    nexttile(i) % rendered structure
    
    imshow(Aggs(inds(i)).image)
    axis('on', 'image')
    set(gca, 'xtick',[], 'xticklabel', [], 'ytick',[], 'yticklabel',[])

    title(titex{i}, 'FontSize', 22, 'interpreter','latex')
end


end

function dist = DISTPP(pp_d, n_bin)
% Get a histogram of primary particle size distributions within an...
%   ...aggregate.
% pp_d: primary particle diameters
% n_bin: number of bins for the histogram
% dist: a structure for primary particle distribution data

n_pp = length(pp_d); % number of primaries
dn_dlogdpp = zeros(n_bin, 1); % discretized size distribution array
fn = zeros(n_bin, 1); % number frequency in each bin
n = zeros(n_bin, 1); % particle counts in each bin
pp_i = zeros(n_pp,1); % initialize placeholder for primaries' bin index
del_dpp = [min(pp_d), max(pp_d)]; % size range

% set the bin locations
r_d = (del_dpp(2) / del_dpp(1))^(1 / n_bin);
d_bin = del_dpp(1) * ones(n_bin + 1, 1);
for i = 1 : n_bin
    d_bin(i+1) = d_bin(i+1) * r_d^(i);
end
d_c = zeros(n_bin,1); % bin centers

for i = 1 : n_bin
    % bin the primary particle sizes 
    if i == 1 % The first bin (only upper limit)
        ii = pp_d < d_bin(i+1);
    elseif i < n_bin
        ii = (pp_d >= d_bin(i)) & (pp_d < d_bin(i+1));
    else % the last bin (all the remaining)
        ii = pp_i == 0;
    end
    
    d_c(i) = sqrt(d_bin(i) * d_bin(i+1)); % get the bin centers
    pp_i(ii) = i; % assign the bin indices to the primaries
    
    n(i) = nnz(i); % count the number of particles in each bin
    fn(i) = n(i) / n_pp; % get the frequency
    dn_dlogdpp(i) = nnz(ii) / log(d_bin(i+1) / d_bin(i)); % calculate the size distribution
end

% assign distribution values obtained to the output structure
dist.n = n;
dist.fn = fn;
dist.dn_dlogdpp = dn_dlogdpp;
dist.pp_i = pp_i;
dist.d_bin = d_bin;
dist.d_c = d_c;

end