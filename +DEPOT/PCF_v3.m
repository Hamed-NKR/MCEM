function [g, r] = PCF_v3(pp, n_r, c_o, c_i, n_g, opts)
% "PCF" calculates the Pair Correlation Function (PCF) of a fractal...
%   ...aggregate based on descretizing primary particles and counting...
%   ...the number of grid points in randomly origined radially swept...
%   ...volume elements. 
% ----------------------------------------------------------------------- %
% 
% Inputs:
%   pp: Primary particle information matrix
%   n_r: Number of radial increments determining the resolution of PCF.
%   c_o: A coefficient determining the maximum extent of computational...
%       ...domain (this is multiplied by the max pairwise primary...
%       ...particle distance).
%   c_i: A coefficient determining the starting point of radial...
%       ...calculation vector (multipled by min primary particle size).
%   n_g: Number of grid points within the smallest primary particle
%   opts: Options for PCF calculations
% ----------------------------------------------------------------------- %
% 
% Outputs:
%   g: Radial vector of averaged (over primary particles) PCF within and...
%       ...nearby an aggregate.
%   r: The raidal logarithmic points of calculation for PCF starting...
%       ...from inside a central primary particle.
% ----------------------------------------------------------------------- %

% make the options variable if not inputted
if ~exist('opts', 'var') 
    opts = struct();
end

% initialize the visibility variable
if (~isfield(opts, 'vis')) || isempty(opts.vis)
    opts.vis = 'on'; % default to plot the results
end

% initialize the textbar display variable
if (~isfield(opts, 'tbar')) || isempty(opts.tbar)
    opts.tbar = 'on'; % default to print the calculation progress
end

% initialize the logging origin variable
if (~isfield(opts, 'orig')) || isempty(opts.orig)
    opts.orig = 'rand'; % default to start from the center of primaries
end

% initialize figure 
if strcmp(opts.vis, 'on') || strcmp(opts.vis, 'ON') || strcmp(opts.vis, 'On')
    figure;
    h = gcf;
    h.Position = [0, 0, 700, 700];
    set(h, 'color', 'white');
end

% initialize resolution parameter
if ~exist('n_r', 'var') || isempty(n_r); n_r = 100; end

% initialize extension parameter
if ~exist('c_o', 'var') || isempty(c_o); c_o = 2; end

% initialize extension parameter
if ~exist('c_i', 'var') || isempty(c_i); c_i = 0.1; end

% initialize grid resolution
if ~exist('n_g', 'var') || isempty(n_g); n_g = 50; end

n_pp = size(pp,1);
inds_pp = nchoosek(1:n_pp, 2);

% maximum pair-wise primary particle distance
r_o = c_o * max(sqrt(sum((pp(inds_pp(:,1),3:5) - pp(inds_pp(:,2),3:5)).^2, 2)) +...
    (pp(inds_pp(:,1),2) + pp(inds_pp(:,2),2)) / 2) / 2;
r_i = c_i * min(pp(:,2));

% initialize pair correlation function vector
r = zeros(n_r + 1,1);
g = zeros(n_r + 1,1);

rr = (r_o / r_i)^(1 / n_r); % radial increment factor
r0 = r_i * ones(n_r + 1, 1); % initialize auxiliary vector for radial incrementation

g0 = zeros(1,n_pp); % placeholder for temporary PCF values from each...
    % ...individual primary particle

% make a grid set within the primaries
r_ppdis = cell(n_pp, 1); % placeholder to store grid points
ind_dmin = find(pp(:,2) == min(pp(:,2)),1);
[r_ppdis{ind_dmin}, c_ppdis] = PAR.MCDISCRETIZEPP_v2(pp(ind_dmin,2),...
    pp(ind_dmin,3:5), n_g); % discretize smallest primary particle

% discretize the rest of primaries based on the concentration of smallest one
for i = 1 : n_pp
    if i~= ind_dmin
        n_ppdis = round(c_ppdis * pi * pp(i,2)^3 / 6);
        r_ppdis{i} = PAR.MCDISCRETIZEPP_v2(pp(i,2), pp(i,3:5), n_ppdis);
    end
end
r_ppdis = cat(1, r_ppdis{:});

% Initialize textbar
if strcmp(opts.tbar, 'on') || strcmp(opts.tbar, 'ON') || strcmp(opts.tbar, 'On')
    fprintf('Volume sweeping started...')
    disp(' ')
    UTILS.TEXTBAR([0, n_r + 1]);
end

for i = 1 : n_r + 1
    % make radial incrementation point
    if i > 1
        r0(i) = r0(i) * rr^(i-1);
        r(i) = sqrt(r0(i-1) * r0(i)); 
    end   
    
    % PCF equals number of grid points within the swept volume over that volume
    for j = 1 : n_pp
        % define the origin of radial volume sweeping
        if strcmp(opts.orig, 'cntr') || strcmp(opts.orig, 'CNTR') || strcmp(opts.orig, 'Cntr')
            rc = pp(j,3:5); % origin to be the center of each primary particle
        elseif strcmp(opts.orig, 'rand') || strcmp(opts.orig, 'RAND') || strcmp(opts.orig, 'Rand')
            % origin to be a random point within each primary particle
            rc0 = rand(1,3); % rho, theta and phi values, respectively, in a cylindrical coordinate
            rc = zeros(1,3); % initialize the origin coordinates
            rc(1) = (pp(j,2) / 2) * rc0(1) .* cos(2 * pi * rc0(2)) .* sin(pi * rc0(3)) + pp(j,3); % x = rho * cos(theta) * sin(phi) 
            rc(2) = (pp(j,2) / 2) * rc0(1) .* sin(2 * pi * rc0(2)) .* sin(pi * rc0(3)) + pp(j,4); % y = rho * sin(theta) * sin(phi) 
            rc(3) = (pp(j,2) / 2) * rc0(1) .* cos(pi * rc0(3)) + pp(j,5); % z = rho * cos(phi)
        end
        
        r_pair = sqrt(sum((r_ppdis - rc).^2, 2));
        
        if i > 1
            g0(j) = nnz((r_pair >= r0(i-1)) & (r_pair < r0(i))) /...
                (4 * pi * r(i)^2 * (r0(i) - r0(i-1))) / c_ppdis;
        else
            g0(j) = nnz(r_pair < r0(i)) /...
                (4 * pi * r0(i)^3) / (3 * c_ppdis);
        end            
    end
    
    g(i) = mean(g0);
    
    if strcmp(opts.tbar, 'on') || strcmp(opts.tbar, 'ON') || strcmp(opts.tbar, 'On')
        UTILS.TEXTBAR([i, n_r + 1]); % Update textbar
    end
end

g = g(g~=0);
r = r(g~=0);
    
% plot PCF vs. normalized radial distance averaged over different primary particles
r_norm = geomean(pp(:,2)) / 2;
r = r / r_norm;

if strcmp(opts.vis, 'on') || strcmp(opts.vis, 'ON') || strcmp(opts.vis, 'On')
    plot(r, g);
    hold on

    box on
    set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 18,...
        'TickLength', [0.02 0.02], 'XScale', 'log', 'YScale', 'log')    
    xlabel('$\overline{r}$ (-)', 'interpreter', 'latex', 'FontSize', 20)
    ylabel('$\overline{g}$($\overline{r}$) (-)', 'interpreter', 'latex', 'FontSize', 20)
end

end

