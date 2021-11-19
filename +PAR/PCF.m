function [g, r, morph, h] = PCF(pars, kk_pars, kk_h)
% "PCF" calculates the Pair Correlation Function (PCF) for single...
%   ...aggregates. The relations used are based on Heinson et al. ...
%   ...(2012)'s work.
% ----------------------------------------------------------------------- %
% 
% Inputs:
%   pars: Particle information structure/class
%   kk_pars: Indices of the aggregates to be analyzed
%   kk_h: ~ to be plotted
% ----------------------------------------------------------------------- %
% 
% Outputs:
%   g: Pair correlation function with respect to the center of mass...
%       ...(COM) within the aggregate domain
%   r: A discuntinuous set of distances from com within the aggregate
%   morph: A data structure containing the fratcal properties of the...
%       ...individual aggregates
%   h: Output figure handle
% ----------------------------------------------------------------------- %

% Total number of aggregates
if isa(pars, 'AGG')
    n_tot = length(pars);
else
    n_tot = length(pars.n);
end

if ~exist('kk_pars', 'var') || isempty(kk_pars)
    kk_pars = 1 : n_tot;
end

if ~exist('kk_h', 'var')
    kk_h = [];
end

% Initializing the figure handle
if isempty(k_hh)
    h = [];
elseif length(kk_h) > 24
    error('Invalid number of output images! (should be <= 24)')
else
    hold off
    h = gcf;

    % Clearing previous data (for animations)
    all_axes_in_figure = findall(h, 'type', 'axes');
    n_ax = numel(all_axes_in_figure);
    for i = 1 : n_ax
        cla(all_axes_in_figure(i))
    end

    figure(h);

    % Setting figure position and background
    h.Position = [0, 0, 2000, 892.1];
    set(h, 'color', 'white');
    
    tiledlayout('flow')
end

% Compiling the aggregates' properties
if isa(pars, 'AGG')
    pp = cat(1, pars(kk_pars).pp);
    n = cat(1, pars(kk_pars).n);
    dmax = cat(1, pars(kk_pars).dmax);
else
    pp = cat(1, pars.pp(kk_pars));
    n = cat(1, pars.n(kk_pars));
    dmax = cat(1, pars.dmax(kk_pars));
end

n_agg = length(kk_pars); % Number of aggregates to be analyzed

r = cell(n_agg, 1); % Initializing the discrete radial distance from COM
g = cell(n_agg, 1); % The corresponding discretized PCF

% Center of mass of each aggregate
if ~exist('pars.r', 'var') || isempty(pars.r)
    r_com = PAR.COM(pp, n); 
else
    if isa(pars, 'AGG')
        r_com = cat(1, pars(kk_pars).r);
    else
        r_com = cat(1, pars.r(kk_pars));
    end
end
r_com = repelem(r_com, n, 1);
r_com = mat2cell(r_com, n);

% Mean Primary particle size within each aggregate
if ~exist('pars.dpp', 'var') || isempty(pars.dpp)
    dpp = PAR.MEANPP(pp); 
else
    if isa(pars, 'AGG')
        dpp = cat(1, pars(kk_pars).dpp);
    else
        dpp = cat(1, pars.dpp(kk_pars));
    end
end
dpp = dpp(:,1); % Just the mean, removing std data

% Initializing the structure containing fractal properties
morph = struct('df', zeros(n_agg,1), 'phi', zeros(n_agg,1),...
    'gamma', zeros(n_agg,1), 'zeta', zeros(n_agg,1));

for i = 1 : n_agg
    
    r{i} = 0 : dmax(i) / (10 * n(i)) : dmax(i);
    dr_pp = sqrt(sum((pp{i}(:, 3:5) - r_com{i}).^2, 2)) - pp{i}(:, 2) / 2;
        % Nearest distance of each primary particle from COM
    g{i} = zeros(1, length(r{i}));
    
    for j = 2 : numel(r{i})
        g{i}(j) = nnz(dr_pp < r{i}(j)) / (4 * pi * r{i}(j)^2 *...
            (r{i}(j) - r{i}(j))); % Number of primaries closer than a...
            % ...certain distance over the swept volume
    end
    
%     [g{i}, iu] = unique(g{i}, 'stable'); % Removing repeated elements
%     r{i} = r{i}(iu);
    
%     rdata = r{i}; % Storing r to be used in the fitting function
    pcfit = @(m, rdata) m(1) * (rdata^m(2)) * exp(m(3) * rdata^m(4));
        % The theoretical PCF function to be fitted (see Heinson...
            % ...et al. (2012))
    dr0 = [0, dmax(i)]; % Range of fit
    m = lsqcurvefit(pcfit, dr0, r{i}, g{i}); % Least-square nonlinear...
        % ...fitter for PCF 
    
    % Fractal properties
    morph.df(i) = m(2) + 3; % Fractal dimension
    morph.gamma(i) = m(4); % Stretching exponent
    morph.zeta(i) = (-m(2))^(-1/m(4)); % Stretching prefactor
    morph.df(i) = m(1) * 4 * pi * (dpp(i)^(m(2) + 3)) / (m(2) + 3);
        % Packing factor
    
    if ismember(kk_pars(i), kk_h)
        nexttile
        
        % PCF vs. radial distance from COM data
        scatter(r{i}, g{i}, 25, [0.8500 0.3250 0.0980], 'filled');
        hold on
        
        % Least square fit
        r_discrete = linspace(r{i}(1), r{i}(end), 1000);
        plot(xfit, pcfit(m, r_discrete), 'Color', [0 0.4470 0.7410],...
            'LineWidth', 2.5);
        
        % Plot specs
        title('Aggregate id: ' + num2str(kk_pars(i), '%d'),...
            'FontName', 'Times New Roman', 'FontWeight', 'bold');
        xlabel('r (m)', 'FontName', 'Times New Roman',...
            'FontWeight', 'bold')
        ylabel('g(r) (m^-3)', 'FontName', 'Times New Roman',...
            'FontWeight', 'bold')
        set(gca, 'FontName', 'Times New Roman')
        legend(['PCF data'; 'LS fit'], 'Location', 'northwest',...
            'FontName', 'Times New Roman')
        annotation('textbox', 'String', 'd_f = ' +...
            num2str(morph.df(i), '%.2f') + '\phi = ' +...
            num2str(morph.phi(i), '%.2f') + '\gamma = ' +...
            num2str(morph.gamma(i), '%.2f') + '\zeta = ' +...
            num2str(morph.zeta(i), '%.2f'), 'FitBoxToText', 'on',...
            'FontName', 'Times New Roman'); % Fractal properties annotated
        axis padded
    end
    
end


if nargout < 4
    clear h; % Deleting figure handle if not requested as an output
end

end

