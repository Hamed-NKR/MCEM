function dist = DISTPP(pp_d, n_bin)
% Get a histogram of primary particle size distributions within an...
%   ...aggregate.
% pp_d: primary particle diameters
% n_bin: number of bins for the histogram
% dist: a structure for primary particle distribution data

n_pp = length(pp_d); % number of primaries
dn_dlogdpp = zeros(n_bin, 1); % discretized size distribution array
fn = zeros(n_bin, 1); % number frequency in each bin
n = zeros(n_bin, 1); % particle counts in each bin
pp_i = zeros(n_pp,1); % initialize placeholder for primaries' bin index
del_dpp = [min(pp_d), max(pp_d)]; % size range

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
dist.pp_i = pp_i;
dist.d_bin = d_bin;
dist.d_c = d_c;

end
