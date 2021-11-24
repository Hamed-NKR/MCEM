function [pa_avg, h] = PROJECTION(pars, rsl_samp, rsl_avg, kk_pars, kk_h)
% "PROJECTION" computes the orientation averaged projected area of the...
%   ...agggregates based on Monte Carlo (MC) approach.
% ----------------------------------------------------------------------- %
% 
% Inputs:
%   pars: Particle information structure/class
%   rsl_samp: Sampling resolution (number of spatial points used by MC)...
%       ...at each view angle for each aggregate
%   rsl_avg: Number of random orientations used for averaging the area...
%       ...for each aggregate in MC framework
%   kk_pars: Indices of the aggregates to be analyzed
%   kk_h: ~ to be plotted
% ----------------------------------------------------------------------- %
%
% Outputs:
%   pa_avg: Average projected area array for each aggregate over diffeerent
%       ...angles
%   h: Output figure handle (binary response plots of the samples points...
%       ...being inside(1)/outside(0) the aggregate layed over the...
%       ...original aggregate outline)
% ----------------------------------------------------------------------- %

% Total number of aggregates
if isa(pars, 'AGG')
    n_tot = length(pars);
else
    n_tot = length(pars.n);
end

% Determining the sampling resolution if missing
if ~exist('rsl_samp', 'var') || isempty(rsl_samp)
    rsl_samp = 1e6;
elseif rsl_samp < 1e2
    error('Spatial resolution too low! (better be >= 100)')
end

% Determining the averaging resolution if missing
if ~exist('rsl_avg', 'var') || isempty(rsl_avg)
    rsl_avg = 1e2;
elseif rsl_avg < 10
    error('Angular resolution too low! (better be >= 10)')
end
jj = round(1 + (rsl_avg - 1) .* (0 : 1/3 : 1)'); % Indices of angles to...
% ...be plotted

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
elseif length(kk_h) > 6
    error('Out of range number of aggregates requested for plotting! (should be <= 6)')
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
    
    tiledlayout(length(kk_h),4);
end

% Compiling the aggregates' properties
if isa(pars, 'AGG')
    pp0 = cat(1, pars(kk_pars).pp);
    n = cat(1, pars(kk_pars).n);
else
    pp0 = cat(1, pars.pp(kk_pars));
    n = cat(1, pars.n(kk_pars));
end

n_agg = length(kk_pars); % Number of aggregates to be analyzed

pa = zeros(n_agg, length(rsl_avg)); % Initializing the projected area set

for i = 1 : n_agg
    
    angs = 2 * pi * rand(rsl_avg, 3); % A uniform random set of three...
        % ...intrinsic Euler angles for rotating the aggregtes to...
        % ...average the projected area
    
    for j = 1 : rsl_avg
        pp = cell2mat(PAR.ROTATE(pp0(i), n(i), angs(j,:))); % Rotating...
            % ...the aggregate in random direction and converting to matrix
        
        x_rng = [min(pp(:,3)), max(pp(:,3))]; % Extension range in x dir.
        y_rng = [min(pp(:,4)), max(pp(:,4))]; % ~ y dir.
        
        pp_temp = repmat(pp, rsl_samp, 1);
        
        mcpoints = [x_rng(1), y_rng(1)] + [x_rng(2) - x_rng(1),...
            y_rng(2) - y_rng(1)] .* rand(rsl_samp, 2); % Locating the...
                % ...MC uniformly random points 
        
        mcstat = sqrt(sum((pp_temp(:, 3:4) -...
            repelem(mcpoints, n(i), 1)).^2, 2)) <= pp_temp(:, 2);
                % Check if the MC points fall within primaries
        mcstat = mat2cell(mcstat, n(i) * ones(rsl_samp,1));
        
        % Determine if the MC points are inside the aggregates
        for k = 1 : rsl_samp
            mcstat{k} = any(mcstat{k});
        end
        mcstat = cell2mat(mcstat);
        
        pa(i,j) = nnz(mcstat) / length(mcpoints); % Fraction of points...
            % ...falling inside the domain
        
        % Plotting the results
        if ismember(kk_pars(i), kk_h) && ismember(j, jj) % Plotting only...
                % ...for 4 angles
            nexttile

            % The original outline of aggregate (i)
            viscircles([pp(:,3), pp(:,4)], pp(:,2) ./ 2,...
                'EnhanceVisibility', false, 'Color', [0 0.4470 0.7410],...
                'LineWidth', 1.5); 
            hold on
                        
            % MC points
            scatter(mcpoints(mcstat,1), mcpoints(mcstat,2), 1,...
                [0.8500 0.3250 0.0980], 'filled');
            % Plot specs
            title("Agg id: " + num2str(kk_pars(i), "%d"),...
                'FontName', 'Times New Roman', 'FontWeight', 'bold');
            xlabel('x (m)', 'FontName', 'Times New Roman',...
                'FontWeight', 'bold')
            ylabel('y(m)', 'FontName', 'Times New Roman',...
                'FontWeight', 'bold')
            xlim(x_rng)
            ylim(y_rng)
            set(gca, 'FontName', 'Times New Roman')
            if find(kk_h == kk_pars(i)) == 1 &&...
                    find(jj == j) == 1 % Only one legend for the whole plot
                lgd = legend({'Agg outline', 'MC points'},...
                'FontName', 'Times New Roman');
                lgd.Layout.Tile = 'south';
            end
            axis padded
        end
    end
end        

pa_avg = mean(pa,2); % Averaging over orientations

if nargout < 2
    clear h; % Deleting figure handle if not requested as an output
end

end

