function [pa_ens, h] = EVALPROJ(pars, kk_pars, n_samp, rsl_samp,...
    rsl_avg)
% "EVALPROJ" examines the performance of the Monte Carlo (MC) algorithm...
%   ...embedded in MCEM to compute aggregate projected area.
% ----------------------------------------------------------------------- %
% 
% Inputs:
%   pars: Particle information structure/class
%   kk_pars: Indices of the aggregates to be analyzed
%   n_samp: Number of tests taken for each aggregate
%   rsl_samp: Sampling resolution (number of spatial points used by MC)...
%       ...at each view angle for each aggregate
%   rsl_avg: Number of random orientations used for averaging the area...
%       ...for each aggregate in MC framework
% ----------------------------------------------------------------------- %
%
% Outputs:
%   pa_avg: Average projected area array for each aggregate over diffeerent
%       ...angles
%   af_avg: Fraction of orientation averaged projected area to the...
%       ...surrounding rectangle
%   h: Output figure handle (statistic results of projected area)
% ----------------------------------------------------------------------- %

% Total number of aggregates
if isa(pars, 'AGG')
    n_tot = length(pars);
else
    n_tot = length(pars.n);
end

% Assigning the aggregates to be analyzed if not given
if ~exist('kk_pars', 'var') || isempty(kk_pars)
    kk_pars = 1 : n_tot;
end

% Assigning the number of MC tests to be analyzed if not given
if ~exist('n_samp', 'var') || isempty(n_samp)
    n_samp = 5;
end

% Determining the sampling resolution if missing
if ~exist('rsl_samp', 'var') || isempty(rsl_samp)
    rsl_samp = 1e4;
elseif rsl_samp < 1e2
    error('Spatial resolution too low! (better be >= 100)')
end

% Determining the averaging resolution if missing
if ~exist('rsl_avg', 'var') || isempty(rsl_avg)
    rsl_avg = 20;
elseif rsl_avg < 10
    error('Angular resolution too low! (better be >= 10)')
end

pa_ens = zeros(length(kk_pars), n_samp, 2); % Projected area database
for i = 1 : n_samp
    fprintf('Test %d:', i);
    disp(newline)

    [pa_ens(:,i,1), pa_ens(:,i,2)] = PAR.PROJECTION(pars, kk_pars,...
        rsl_samp, rsl_avg); % Storing projected area and area ratio for...
            % ...each test
end

figure;
h = gcf;
h.Position = [0, 0, 700 * length(kk_pars) / 4, 700]; % Position and size
set(h, 'color', 'white'); % Background color

xlabs = cell(length(kk_pars), 1);
for i = 1 : length(kk_pars)
     xlabs{i} = num2str(kk_pars(i));
end

boxplot((100 * pa_ens(:,:,2))', 'Labels', xlabs, 'Whisker',  1)
set(gca, 'FontName', 'SansSerif', 'FontSize', 12)
title('Variations of projected area', 'FontName', 'SansSerif',...
    'FontWeight', 'bold', 'FontSize', 18)
subtitle(' ','FontName', 'SansSerif', 'FontSize', 8)
xlabel('Agg id', 'FontName', 'SansSerif', 'FontWeight', 'bold',...
    'FontSize', 14)
ylabel('AF (%)','FontName', 'SansSerif', 'FontWeight', 'bold',...
    'FontSize', 14)

if nargout < 2
    clear h; % Deleting figure handle if not requested as an output
end

end

