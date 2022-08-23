function [h, dist] = DISTPP(pp_d, opts)
% "DISTPP" plots the size distribution of a population of primary
%   ...particles.
% ----------------------------------------------------------------------- %
% 
% Inputs:
%   pp_d: Primary particle size array
%   params: a structure containing parameters to be plotted.
% ----------------------------------------------------------------------- %
%
% Outputs:
%   h: output figure handle
%   dist: a structure containing distribution data
% ----------------------------------------------------------------------- %

% initialize plotting options variable
if ~exist('opts', 'var') 
    opts = struct();
end

% set number of bins for data fitting
if ~isfield(opts, 'n_bin')
    opts.n_bin = [];
end
n_bin = opts.n_bin;
if isempty(n_bin)
    n_bin = 2^4; % default number of bins
end

% set the primary particle size range
if ~isfield(opts, 'dpp')
    opts.dpp = [];
end
del_dpp = opts.dpp;
if isempty(del_dpp)
    del_dpp = [min(pp_d), max(pp_d)]; % default size range
end

n_pp = length(pp_d); % total number of primaries to be studied

% initialize distribution variables
dn_dlogdpp = zeros(n_bin, 1); % discretized size distribution array
fn = zeros(n_bin, 1); % number frequency in each bin
n = zeros(n_bin, 1); % particle counts in each bin

pp_i = zeros(n_pp,1); % initialize placeholder for primaries' bin index

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

% initialize figure properties
figure
h = gcf;
if ~all(h.Position == [0, 0, 700, 700])
    h.Position = [0, 0, 700, 700];
end
set(h, 'color', 'white')

% plot discretized size distribution
x1 = 1e9 * [d_bin(1); repelem(d_bin(2 : end - 1), 2); d_bin(end)];
y1 = repelem(dn_dlogdpp, 2);
plot(x1, y1, 'Color', [0.1 0.1 0.1], 'LineWidth', 1)

box on
set(gca, 'FontName', 'Calibri Light', 'FontSize', 12, 'TickLength', [0.01 0.01])
set(gca, 'XScale', 'log')
xlabel('{\it{d_{pp}}} [nm]', 'FontName', 'Calibri', 'FontSize', 16)
xlim(1e9 * [del_dpp(1), del_dpp(2)])
ylabel('{d\it{n_{pp}}\rm / d(log\it{d_{pp}}}) [{m}^{-1}]',...
    'FontName', 'Calibri', 'FontSize', 16)
% ylim([0, ceil(max(ylim) / 0.1) * 0.1])


end