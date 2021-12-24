function parsfilter = FILTERPARS(parsdata, kk, opts)
% "FILTERPARS" filters the aggregate data in order to get a more...
    % ...realistic, i.e. random, representation of the aggregate...
    % ...population's fractality.
% ----------------------------------------------------------------------- %
%
% Inputs:
%   parsdata: Stored particle information over time
%   kk: Indices of stored data to be analyzed
%   opts: Input options(type of filter applied)
% ----------------------------------------------------------------------- %
%
% Output:
%   parsfilter: Filtered particle data
% ----------------------------------------------------------------------- %

% Setting options
if ~(exist('opts', 'var') && isfield(opts, 'method'))
    opts = struct('method', 'mean'); % Setting default method to mean
end

if ~ismember(opts.method, {'mean', 'median'})
    error('Invalid particle data filtering method!')
end

% Extracting desired aggregates data
dg_dpp0 = cat(1, parsdata.dg_dpp{kk}); % Ratio of gyration diameter to...
    % ...diameter of primaries
npp0 = cat(1, parsdata.npp{kk}); % Number of primaries

npp = unique(npp0, 'sorted'); % Sorting aggregates based on the number...
    % ...of primaries
n_unq = length(npp); % Number of discrete bins
dg_dpp = cell(n_unq,1); % Cell array to store binned size ratio
dg_dpp_avg = zeros(n_unq,1); % Array to store average values over bins

for i = 1 : n_unq
    dg_dpp{i} = dg_dpp0(npp0 == npp(i));
    switch opts.method
        case 'mean'
            dg_dpp_avg(i) = mean(dg_dpp{i});
        case 'median'
            dg_dpp_avg(i) = median(dg_dpp{i});
    end
end

parsfilter = struct('npp', npp, 'dg_dpp', dg_dpp_avg);

end

