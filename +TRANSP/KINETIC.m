function [beta, tau, z] = KINETIC(n, t, opts)
% "KINETIC" calculates the properties related to kinetics of aggrgeation.
% ----------------------------------------------------------------------- %
%
% Inputs:
%   n: Temporal variations of total number of independent particles
%   t: Times associated with n series
%   opts: plotting options
% ----------------------------------------------------------------------- %
%
% Output:
%   beta: Collision frequency (a.k.a coagulation kernel)
%   tau: Characteristic time of aggregation
%   z: Kinetic exponent
% ----------------------------------------------------------------------- %

% Set plotting defaults
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

ninv =  1./n - 1/n(1); % Data to be fitted for characteristic time and...
    % ...kinetic exponent

n_dat0 = round(length(t) / 20); % timespot to start with for fitting 

% Fit a linear regression model to n inverse time series
fit1 = fitlm(table(log(t(n_dat0 : end)), log(ninv(n_dat0 : end))), 'linear');
z = fit1.Coefficients.Estimate(2);
tau = exp(-fit1.Coefficients.Estimate(1) / z);

% Get the collision kernel
beta = (2 * z / tau^z) * t.^(z - 1);

if ismember(opts_visual, {'ON', 'On', 'on'})
    % initialize figure properties
    figure
    h = gcf;
    if ~all(h.Position == [0, 0, 1200, 600])
        h.Position = [0, 0, 1200, 600];
    end
    set(h, 'color', 'white')
    
    % set the layout
    tt = tiledlayout(1, 2);
    tt.TileSpacing = 'compact';
    tt.Padding = 'compact';
    
    % time history of agg numbers
    tt1 = nexttile;
    p11 = plot(t, ninv, 'Color', [0.1 0.1 0.1], 'LineWidth', 1);
    hold on
    
    % The fit line based on LRM results
    t_fit = (min(t) : range(t) / (10 * (numel(t) - 1)) : max(t))';
    ninv_fit = (t_fit / tau) .^ z;    
    p12 = plot(t_fit, ninv_fit, 'Color', [0.4660 0.6740 0.1880],...
        'LineWidth', 2, 'LineStyle', '-.');
    
    box on
    set(gca, 'XScale', 'log')
    set(gca, 'YScale', 'log')
    set(gca, 'FontName', 'SansSerif', 'FontSize', 12, 'TickLength', [0.02 0.02])
    xlabel('t (s)', 'FontName', 'SansSerif', 'FontWeight', 'bold',...
        'FontSize', 14)
    ylabel('1/n - 1/n_0 (-)', 'FontName', 'SansSerif', 'FontWeight', 'bold', 'FontSize', 14)
    legend([p11, p12], {'Real-time', 'Linear regression fit'},...
        'Location', 'northwest', 'FontName', 'SansSerif', 'FontSize', 12);
    title(tt1, 'Charactertic kinetic time and exponent', 'FontName', 'SansSerif',...
        'FontWeight', 'bold', 'FontSize', 16)
    
    % Collision kernel time history
    tt2 = nexttile;
    beta_hr = (2 * z / tau^z) * t_fit.^(z - 1); % Draw in high resolution
    plot(t_fit, beta_hr, 'Color', [0.1 0.1 0.1], 'LineWidth', 1)
    
    box on
    set(gca, 'XScale', 'log')
    set(gca, 'FontName', 'SansSerif', 'FontSize', 12, 'TickLength', [0.02 0.02])
    xlabel('t (s)', 'FontName', 'SansSerif', 'FontWeight', 'bold',...
        'FontSize', 14)
    ylabel('beta (s^-^1)', 'FontName', 'SansSerif', 'FontWeight', 'bold', 'FontSize', 14)
    title(tt2, 'Collision kernel', 'FontName', 'SansSerif',...
        'FontWeight', 'bold', 'FontSize', 16)    
end

