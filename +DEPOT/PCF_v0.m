function [g, r, morph, h] = PCF_v0(pars, kk_pars, kk_h)
% "PCF" calculates the Pair Correlation Function (PCF) and fractal...
%   ...properties for individual aggregates. The relations used are...
%   ...based on Heinson et al. (2012)'s work.
% ----------------------------------------------------------------------- %
% 
% Inputs:
%   pars: Particle information structure/class
%   kk_pars: Indices of the aggregates to be analyzed
%   kk_h: ~ to be plotted
% ----------------------------------------------------------------------- %
% 
% Outputs:
%   g: Pair correlation function within the aggregate domain...
%       ...with respect to its center of mass (COM) 
%   r: A discrete set of distances from COM for PCF calculations
%   morph: A data structure containing the fratcal properties of the...
%       ...individual aggregates
%   h: Output figure (g vs. r for different aggregates) handle
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

% Initializing the aggregate indices for plotting if not defined
if ~exist('kk_h', 'var')
    kk_h = [];
end

% Initializing the figure handle
if isempty(kk_h)
    h = [];
elseif length(kk_h) > 24
    error('Out of range number of output images! (should be <= 24)')
else
    hold off
    h = gcf;

    % Clearing previous data
    all_axes_in_figure = findall(h, 'type', 'axes');
    n_ax = numel(all_axes_in_figure);
    for i = 1 : n_ax
        cla(all_axes_in_figure(i))
    end

    figure(h);

    % Setting figure position and background
%     h.Position = [0, 0, 2000, 892.1];
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

r = cell(n_agg, 1); % Initializing the discrete radial distance
g = cell(n_agg, 1); % ~ the corresponding discretized PCF

% Center (of mass) of each aggregate
if ~exist('pars.r', 'var') || isempty(pars.r)
    r_c = PAR.COM(pp, n); 
else
    if isa(pars, 'AGG')
        r_c = cat(1, pars(kk_pars).r);
    else
        r_c = cat(1, pars.r(kk_pars));
    end
end
r_c = repelem(r_c, n, 1);
r_c = mat2cell(r_c, n);

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
    dr_pp0 = sqrt(sum((pp{i}(:, 3:5) - r_c{i}).^2, 2)); % Central...
        % ...distance of each primary particle from COM
    j0 = dr_pp0 == min(dr_pp0); % Identifying the nearest primary to COM
    dr_pp = sqrt(sum((pp{i}(:, 3:5) -...
        repmat(pp{i}(j0, 3:5), n(i), 1)).^2, 2)) - pp{i}(:, 2) / 2;
        % Nearest distance of each primary from the center of the...
            % ...nearest primary to COM!!
        
    g{i} = zeros(1, length(r{i}));
%     g{i}(1) = Inf; % Zero in the denominator at the center while the...
%         % ...numerator is 1
    
    for j = 2 : numel(r{i})
        g{i}(j) = nnz(dr_pp < r{i}(j)) / (4 * pi * r{i}(j)^2 *...
            (r{i}(j) - r{i}(j-1))); % Number of primaries closer than...
            % ...the marched distance over the swept volume
    end
    
    % Removing the central point
    r{i} = r{i}(2:end);
    g{i} = g{i}(2:end);
    
%     [g{i}, iu] = unique(g{i}, 'stable'); % Removing repeated elements
%     r{i} = r{i}(iu);
    
%     rdata = r{i}; % Storing r to be used in the fitting function
    pcfit = @(m, rr) m(1) * (rr.^m(2)) .* exp(m(3) * rr.^m(4));
        % The theoretical PCF function to be fitted (see Heinson...
            % ...et al. (2012))
    m0 = [30, -5, -5, 10]; % Initial guess
    m = lsqcurvefit(pcfit, m0, r{i}, g{i}); % Least-square nonlinear...
        % ...fitter for PCF 
    
    % Fractal properties
    morph.df(i) = m(2) + 3; % Fractal dimension
    morph.gamma(i) = m(4); % Stretching exponent
    morph.zeta(i) = (-m(3))^(-1/m(4)); % Stretching prefactor
    morph.df(i) = m(1) * 4 * pi * (dpp(i)^(m(2) + 3)) / (m(2) + 3);
        % Packing factor
    
    % Plotting the results
    if ismember(kk_pars(i), kk_h)
        nexttile
        
        % PCF vs. radial distance from the central primary particle
        scatter(r{i}, g{i}, 25, [0.8500 0.3250 0.0980], 'filled');
        hold on
        
        % Least square fit
        r_discrete = linspace(r{i}(1), r{i}(end));
        plot(r_discrete, pcfit(m, r_discrete), 'Color',...
            [0 0.4470 0.7410], 'LineWidth', 2);
        
        % Set axis scales to be logarithmic
        set(gca, 'XScale', 'log')
        set(gca, 'YScale', 'log')
        
        % Plot specs
        title("Agg id: " + num2str(kk_pars(i), "%d"),...
            'FontName', 'Times New Roman', 'FontWeight', 'bold');
        xlabel('r (m)', 'FontName', 'Times New Roman',...
            'FontWeight', 'bold')
        ylabel('g(r) (m^-^3)', 'FontName', 'Times New Roman',...
            'FontWeight', 'bold')
        set(gca, 'FontName', 'Times New Roman')
        if find(kk_h == kk_pars(i)) == 1 % Only one legend for the whole...
            % ...plot
            lgd = legend({'PCF data', 'LS fit'},...
            'FontName', 'Times New Roman');
            lgd.Layout.Tile = 'south';
        end        
%         anotxt =  "n_pp = " + num2str(n(i), '%d') + newline +...
%             "d_pp = " + num2str(dpp(i,1), '%1.1e') + " (m)" + newline +...
%             "d_f = " + num2str(morph.df(i), '%.2f') + newline +...
%             "\phi = " + num2str(morph.phi(i), '%.2f') + newline +...
%             "\gamma = " + num2str(morph.gamma(i), '%.2f') + newline +...
%             "\zeta = " + num2str(morph.zeta(i), '%.2f');
%         annotation('textbox', 'String', anotxt, 'FitBoxToText', 'on',...
%             'FontName', 'Times New Roman');
%             % Fractal properties annotated
        axis padded
    end
    
end


if nargout < 4
    clear h; % Deleting figure handle if not requested as an output
end

end

