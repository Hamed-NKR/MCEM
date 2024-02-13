function [g, r] = PCF_v3(pp, n_r, c_o, c_i, n_g)
% "PCF" calculates the Pair Correlation Function (PCF) of a fractal...
%   ...aggregate based on descretizing primary particles and counting...
%   ...the number of grid points in randomly origined radially swept...
%   ...volume elements. 
% ----------------------------------------------------------------------- %
% 
% Inputs:
%   pp: Primary particle information matrix
%   n_r: Number of radial calculation points (...-1) determining the...
%       ...resolution of PCF.
%   c_o: A coefficient determining the maximum extent of computational...
%       ...domain (this is multiplied by the max pairwise primary...
%       ...particle distance).
%   c_i: A coefficient determining the starting point of radial...
%       ...calculation vector (multipled by min primary particle size).
%   n_g: Number of grid points within the smallest primary particle
% ----------------------------------------------------------------------- %
% 
% Outputs:
%   g: Radial vector of averaged (over primary particles) PCF within and...
%       ...nearby an aggregate.
%   r: The raidal points of calculation for PCF starting from a random...
%       ...within the distances from pp centers for PCF calculations.
% ----------------------------------------------------------------------- %

% initialize figure 
figure;
h = gcf;
h.Position = [0, 0, 700, 700];
set(h, 'color', 'white');

% initialize resolution parameter
if ~exist('n_r', 'var') || isempty(n_r); n_r = 200; end

% initialize extension parameter
if ~exist('c_o', 'var') || isempty(c_o); c_o = 2; end

% initialize extension parameter
if ~exist('c_i', 'var') || isempty(c_i); c_i = 0.2; end

% initialize grid resolution
if ~exist('n_g', 'var') || isempty(n_g); n_g = 200; end

n_pp = size(pp,1);
inds_pp = nchoosek(1:n_pp, 2);

% maximum pair-wise primary particle distance
r_o = c_o * max(sqrt(sum((pp(inds_pp(:,1),3:5) - pp(inds_pp(:,2),3:5)).^2, 2)) +...
    (pp(inds_pp(:,1),2) + pp(inds_pp(:,2),2)) / 2) / 2;
r_i = c_i * min(pp(:,2));

% initialize pair correlation function vector
r = zeros(n_r,1);
g = zeros(n_r,1);

rr = (r_o / r_i)^(1 / n_r); % radial increment factor
r0 = r_i * ones(n_r + 1, 1); % initialize auxiliary vector for radial incrementation

g0 = zeros(1,n_pp); % placeholder for temporary PCF values from each...
    % ...individual primary particle

% make a grid set within the primaries
r_ppdis = cell(n_pp, 1); % placeholder to store grid points
ind_dmin = find(pp(:,2) == min(pp(:,2)),1);
[r_ppdis{ind_dmin}, c_ppdis] = PAR.MCDISCRETIZEPP(pp(ind_dmin,2),...
    pp(ind_dmin,3:5), n_g); % discretize smallest primary particle

% discretize the rest of primaries based on the concentration of smallest one
for i = 1 : n_pp
    if i~= ind_dmin
        n_ppdis = round(c_ppdis * pi * pp(i,2)^3 / 6);
        r_ppdis{i} = PAR.MCDISCRETIZEPP(pp(i,2), pp(i,3:5), n_ppdis);
    end
end
r_ppdis = cat(1, r_ppdis{:});

% Initialize textbar
fprintf('Volume sweeping started...')
disp(' ')
UTILS.TEXTBAR([0, n_r]);

for i = 2 : n_r + 1
    % make radial incrementation point
    r0(i) = r0(i) * rr^(i-1);
    r(i-1) = sqrt(r0(i-1) * r0(i)); 
    
    for j = 1 : n_pp
        % PCF equals number of grid points within the swept volume over that volume
        r_pp = sqrt(sum((r_ppdis - pp(j,3:5)).^2, 2));
        g0(j) = nnz((r_pp >= r0(i-1)) & (r_pp < r0(i))) /...
            (4 * pi * r(i-1)^2 * (r0(i) - r0(i-1))) / c_ppdis;
    end
    
    g(i-1) = mean(g0);
    UTILS.TEXTBAR([i-1, n_r]); % Update textbar
end

g = g(g~=0);
r = r(g~=0);

    
% plot PCF vs. normalized radial distance averaged over different primary particles
r_norm = geomean(pp(:,2)) / 2;
plot(r / r_norm, g);
hold on

box on
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 18,...
    'TickLength', [0.02 0.02], 'XScale', 'log', 'YScale', 'log')    
xlabel('$\overline{r}$ (-)', 'interpreter', 'latex', 'FontSize', 20)
ylabel('$\overline{g}$($\overline{r}$) (-)', 'interpreter', 'latex', 'FontSize', 20)
    
end

