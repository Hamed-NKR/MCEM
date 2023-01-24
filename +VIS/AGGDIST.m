function h = AGGDIST(pp0, mu_da, std_da, nagg2, opts)
% "AGGDIST" displays the prodcures the aggregate population undergoes...
%   ...to get prepared for post-flame mixing.
% ----------------------------------------------------------------------- %
% 
% Inputs:
%   pp0: primary particle info library
%   mu_da: mean diameter of aggregate distribution
%   std_da: standard deviation of aggregate size
%   nagg2: post-sampling total number of aggregates
%   opts: plotting options
% ----------------------------------------------------------------------- %
%
% Outputs:
%   h: output figure handle
% ----------------------------------------------------------------------- %

if ~exist('opts', 'var')
    opts = struct();
end

% set parameters for projected area calculation
if ~isfield(opts, 'n_pnt') || isempty(opts.n_pnt)
    n_pnt = 1e4; % number of Monte Carlo sampling points
else
    n_pnt = opts.n_pnt;
end
if ~isfield(opts, 'n_ang') || isempty(opts.n_ang)
    n_ang = 20; % number of orientations to average over
else
    n_ang = opts.n_ang;
end

pp0 = pp0(:); % Convert the pp data into a 1d cell array if not already

pp0 = pp0(~cellfun('isempty', pp0)); % Remove unused cells if any

pp0 = cat(1, pp0(:)); % Merge pp info from different times

% get number of primaries within aggregates
nagg0 = numel(pp0);
npp0 = zeros(nagg0,1);
dpp0 = zeros(nagg0,1);
for i = 1 : nagg0
    npp0(i) = size(pp0{i}, 1);
    dpp0(i) = geomean(pp0{i}(:,2));
end

pars.pp = pp0;
pars.n = npp0;
da0 = 2 * sqrt(PAR.PROJECTION(pars, [], n_pnt, n_ang) / pi); % get projected area diameter before rescaling

D_TEM = 0.35; % polydispersity exponent
pp1 = pp0;
dpp_emh = ((17.8^(1 / D_TEM) / 100) * (npp0 / 1.1).^...
    (1 / (2 * 1.08))).^(D_TEM / (1 - D_TEM)); % mean pp diameters after rescaling
rpp = 1e-9 * dpp_emh ./ dpp0; % scaling ratios
dpp1 = zeros(nagg0,1);
for i = 1 : nagg0
    pp1{i}(:,2:5) = pp1{i}(:,2:5) * rpp(i); % apply rescaling
    dpp1(i) = geomean(pp1{i}(:,2));
end

% prepare a structure for projected area calculations
pars.pp = pp1;
pars.n = npp0;

da1 = 2 * sqrt(PAR.PROJECTION(pars, [], 1e4, 20) / pi); % get projected area diameter after rescaling

% set aggregate size bins
if ~isfield(opts, 'n_bin') || isempty(opts.n_bin)
    n_bin = 20;
else
    n_bin = opts.n_bin;
end
da_bnd = zeros(2,1);
if ~isfield(opts, 'da_bnd') || isempty(opts.da_bnd)
    [da_bnd(1), da_bnd(2)] = bounds(da1);
else
    da_bnd = opts.da_bnd;
end
r_da = (da_bnd(2) / da_bnd(1))^(1 / n_bin);
da_bin = da_bnd(1) * ones(n_bin + 1, 1); % bin ranges
for i = 1 : n_bin
    da_bin(i+1) = da_bin(i+1) * r_da^(i);
end

% initialize distribution variables
na1 = zeros(n_bin,1); % pre-sampling bin counts
dna1_dlogda = zeros(n_bin,1); % pre-sampling size distribution
nfit = zeros(n_bin,1); % bin counts for the distribution fit
fnfit = zeros(n_bin,1); % number frequency to be fitted
na2 = zeros(n_bin,1); % post-sampling bin counts
dna2_dlogda = zeros(n_bin,1); % post-sampling size distribution

% assign aggregates to the bins
for i = 1 : n_bin
    fnfit(i) = normcdf(log(da_bin(i+1)), log(mu_da), log(std_da)) -...
        normcdf(log(da_bin(i)), log(mu_da), log(std_da)); % the frequency of...
        % ...lognormal fit to the size distribution data
    
    if i == 1 % The first bin (only upper limit)
        ii = da1 < da_bin(i+1);
    elseif i < n_bin
        ii = (da1 >= da_bin(i)) & (da1 < da_bin(i+1)); % middle bins
    else % the last bin (all the remaining)
        ii = da1 >= da_bin(i);
    end
    
    na1(i) = nnz(ii);
    dna1_dlogda(i) = na1(i) / log(da_bin(i+1) / da_bin(i));
end

if ~exist('nagg2', 'var') || isempty(nagg2)
    imax = find(mu_da < da_bin, 1) - 1;
    nagg2 = ceil(na1(imax) / fnfit(imax));
end

% apply random lognormal filter
for i = 1 : n_bin
    nfit(i) = ceil(fnfit(i) * nagg2); % number of particles to be...
        % ...selected in each bin
    na2(i) = nfit(i);
    
%     if nfit(i) <= na1(i) % if enough particles exist in the bin for random sampling
%         na2(i) = nfit(i);
%     
%     elseif na1(i) > 0 % if particle counts in the bin are insufficient
%         na2(i) = ceil(nfit(i) / na1(i)) * na1(i); % duplicate particles
%     end
    
    dna2_dlogda(i) = na2(i) / log(da_bin(i+1) / da_bin(i)); 
end

% initialize figure 
figure;
h = gcf;
h.Position = [0, 0, 1200, 600];
set(h, 'color', 'white');

% initialize layout
tt = tiledlayout(1,2);
tt.TileSpacing = 'loose';
tt.Padding = 'compact';

% display rescaling effect
nexttile

% plot universal correlation
dpp_uc = logspace(0, 4, 1000);
da_uc = 100 * (dpp_uc / 17.8).^(1 / 0.35);
p11 = plot(da_uc, dpp_uc, 'Color', [0.4940 0.1840 0.5560],...
    'LineStyle', '-.', 'LineWidth', 3);
hold on

% plot before scaling population
p12 = scatter(1e9 * da0, 1e9 * dpp0, 10, [0.6, 0.6, 0.6], 'o',...
    'LineWidth', 1);

% plot after scaling population
p13 = scatter(1e9 * da1, 1e9 * dpp1, 10, [0.1, 0.1, 0.1], '^',...
    'LineWidth', 1);

box on
set(gca, 'FontSize', 18, 'TickLength', [0.02 0.02], 'XScale', 'log',...
    'YScale', 'log', 'TickLabelInterpreter','latex')
xlabel('$d_a$ [nm]', 'FontSize', 20, 'interpreter','latex')
xlim([35, 2500])
ylim([10, 100])
ylabel('$\overline{d}_{pp}$ [nm]', 'FontSize', 20, 'interpreter', 'latex')
legend([p12, p13, p11], {'Raw DLCA library', 'Rescaled aggregates', 'Universal correlation'},...
    'Orientation', 'horizontal', 'Location', 'northoutside', 'FontSize', 16,...
    'NumColumns', 2, 'interpreter', 'latex');

% display lognormal sampling outcome
nexttile

x_fit = (1e9 * exp(linspace(log(da_bnd(1)), log(da_bnd(2)), 1000)))';
y_fit = nagg2 * pdf('Normal', log(x_fit), log(1e9 * mu_da), log(std_da));
x_bin = 1e9 * [da_bin(1); repelem(da_bin(2 : end - 1), 2); da_bin(end)];
y_bin1 = repelem(dna1_dlogda, 2);
y_bin2 = repelem(dna2_dlogda, 2);

% plot the target lognormal distribution
p21 = plot(x_fit, y_fit, 'Color', [0.4940 0.1840 0.5560], 'LineWidth', 3,...
    'LineStyle', '-.');
hold on
p22 = plot(x_bin, y_bin1, 'Color', [0.6 0.6 0.6], 'LineWidth', 1.5,...
    'LineStyle', ':');
p23 = plot(x_bin, y_bin2, 'Color', [0.1 0.1 0.1], 'LineWidth', 1.5);

box on
set(gca, 'FontSize', 18, 'TickLength', [0.02 0.02], 'XScale', 'log',...
    'TickLabelInterpreter','latex')
xlabel('$d_a$ [nm]', 'FontSize', 20, 'interpreter','latex')
xlim([min(x_fit), max(x_fit)])
ylim([0, 2200])
ylabel('d$n_a$/dlog($d_a$) [-]', 'FontSize', 20, 'interpreter', 'latex')
legend([p22, p23, p21], {'Unfiltered population', 'Filtered sample', 'Target distribution'},...
    'Orientation', 'horizontal', 'Location', 'northoutside', 'FontSize', 16,...
    'NumColumns', 2, 'interpreter', 'latex');

end

