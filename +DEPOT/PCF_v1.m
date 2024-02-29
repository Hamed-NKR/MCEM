function [g, r, morph, h] = PCF_v1(pp, kk, f_rsl)
% "PCF" calculates the Pair Correlation Function (PCF) and fractal...
%   ...properties for individual aggregates. The relations used are...
%   ...based on Heinson et al. (2012)'s work.
% ----------------------------------------------------------------------- %
% 
% Inputs:
%   pp: A cell array of primary particle information
%   kk: Indices of the aggregates to be plotted
%   f_rsl: spatial resolution factor
% ----------------------------------------------------------------------- %
% 
% Outputs:
%   g: Average pair correlation function within the aggregate domains
%   r: A discrete set of raidal distances from pp centers for PCF calculations
%   morph: A data structure containing the fratcal properties of the...
%       ...individual aggregates
%   h: Output figure (g vs. r for different aggregates) handle
% ----------------------------------------------------------------------- %

% Number of aggregates to be analyzed
n_agg = length(pp);

% Initializing the aggregate indices for plotting if not defined
if ~exist('kk', 'var')
    kk = [];
end

% Initializing the figure handle
if isempty(kk)
    h = [];
elseif length(kk) > 24
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

% resolution parameters
if ~exist('f_rsl', 'var') || isempty(f_rsl); f_rsl = [100, 10]; end

r = cell(n_agg, 1); % Initializing the discrete radial distance
g = cell(n_agg, 1); % ~ the corresponding discretized PCF

% Initializing the structure containing fractal properties
morph = struct('df', zeros(n_agg,1), 'phi', zeros(n_agg,1),...
    'gamma', zeros(n_agg,1), 'zeta', zeros(n_agg,1));

% Compiling the aggregates' properties
n_pp = zeros(n_agg,1); % number of primary particles within the aggregates
r_agg = zeros(n_agg,2); % min and max radial distance within the aggregates
for i = 1 : n_agg
    n_pp(i) = size(pp{i},1);
    r_agg(i,1) = min(pp{i}(:,2)) / 2 / f_rsl(2);
    ind_pp = nchoosek(1 : n_pp, 2);
    
    % maximum pair-wise primary particle distance
    r_agg(i,2) = f_rsl(2) * max(sqrt(sum((pp{i}(ind_pp(:,1),3:5) - pp{i}(ind_pp(:,2),3:5)).^2, 2)) +...
        (pp{i}(ind_pp(:,1),2) + pp{i}(ind_pp(:,2),2)) / 2) / 2; 
        
    r{i} = zeros(f_rsl(1) * n_pp(i) ,1);
    g{i} = zeros(f_rsl(1) * n_pp(i) ,1);
    
    rr = (r_agg(i,2) / r_agg(i,1))^(1 / (f_rsl(1) * n_pp(i)));
    r0 = r_agg(i,1) * ones(f_rsl(1) * n_pp(i) + 1, 1);
    
    g0 = zeros(1,n_pp(i));
    
    for j = 2 : f_rsl(1) * n_pp(i) + 1
        r0(j) = r0(j) * rr^(j-1);
        r{i}(j-1) = sqrt(r0(j-1) * r0(j));
        
        for k = 1 : n_pp(i)
            % min distance of pp edges from the center of primary particle of interest
            r_pp = sqrt(sum((pp{i}(:,3:5) - repmat(pp{i}(k,3:5), n_pp(i), 1)).^2, 2)) - ...
                pp{i}(:,2)/2;            

            % Number of primaries closer than the marched distance over the swept volume
            if r0(j) <= pp{i}(k,2)/2
                g0(k) = 1 / (pi * pp{i}(k,2)^3);
            else
                dn = nnz((r_pp <= repmat(r0(j-1), n_pp(i), 1)) &...
                    (r_pp + pp{i}(:,2)) > repmat(r0(j), n_pp(i), 1)) ;
                g0(k) =  dn/ (4 * pi * r{i}(j-1)^2 * (r0(j) - r0(j-1)));
            end
        end
        
        g{i}(j-1) = mean(g0);
    end
    
    g{i} = g{i}(g{i}~=0);
    r{i} = r{i}(g{i}~=0);
    
%     pcfit = @(m, rf) m(1) * (rf.^m(2)) .* exp(m(3) * rf.^m(4));
%         % The theoretical PCF function to be fitted (see Heinson...
%             % ...et al. (2012))
%     m0 = [30, -5, -5, 10]; % Initial guess
%     m = lsqcurvefit(pcfit, m0, r{i}, g{i}); % Least-square nonlinear...
%         % ...fitter for PCF 
%     
%     % Fractal properties
%     morph.df(i) = m(2) + 3; % Fractal dimension
%     morph.gamma(i) = m(4); % Stretching exponent
%     morph.zeta(i) = (-m(3))^(-1/m(4)); % Stretching prefactor
%     morph.df(i) = m(1) * 4 * pi * (dpp(i)^(m(2) + 3)) / (m(2) + 3);
%         % Packing factor
    
    % Plotting the results
    if ismember(i, kk)
        nexttile
        
        % PCF vs. radial distance from the central primary particle
        plot(1e9 * r{i}, 1e-27 * g{i});
        hold on
        
%         % Least square fit
%         r_discrete = linspace(r{i}(1), r{i}(end));
%         plot(r_discrete, pcfit(m, r_discrete), 'Color',...
%             [0 0.4470 0.7410], 'LineWidth', 2);
        
        % Set axis scales to be logarithmic
        set(gca, 'XScale', 'log')
        set(gca, 'YScale', 'log')
        
        % Plot specs
%         title("Agg id: " + num2str(kk(i), "%d"),...
%             'FontName', 'Times New Roman', 'FontWeight', 'bold');
        xlabel('r (nm)', 'FontName', 'Times New Roman',...
            'FontWeight', 'bold')
        ylabel('g(r) (nm^-^3)', 'FontName', 'Times New Roman',...
            'FontWeight', 'bold')
        set(gca, 'FontName', 'Times New Roman')
%         if find(kk == kk_pars(i)) == 1 % Only one legend for the whole...
%             % ...plot
%             lgd = legend({'PCF data', 'LS fit'},...
%             'FontName', 'Times New Roman');
%             lgd.Layout.Tile = 'south';
%         end        
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

