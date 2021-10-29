function parsdata = SAVEPARS(pars, time, k, params_ud, parsdata)
% "SAVEPARS" records the particle data at different timesteps in a...
%   ...structural format with different fields of properties.
% ----------------------------------------------------------------------- %
% 
% Inputs/Output:
%   pars: Particle information structure/class
%   time: Time array
%   k: Iteration number
%   parsdata: Structure for particle data to be saved in
% ----------------------------------------------------------------------- %

% Initializing the output structure
if ~exist('parsdata', 'var')
    parsdata = struct('ii', [], 't', [], 'ntot', [], 'vf', [],...
        'beta', [], 'dn_dlogdv', [], 'dm_dlogdv', [], 'dg_dpp', [],...
        'npp',[], 'df', [], 'kf', []);
end

% Compiling properties across aggregates
dv = cat(1, pars.dv);
m = cat(1, pars.m);
dg = cat(1, pars.dg);
dpp = cat(1, pars.dpp);
npp = cat(1, pars.n);

parsdata.ii = [parsdata.ii; k]; % Adding iteration index
parsdata.t = [parsdata.t; time(k)]; % Time

if isa(pars, 'AGG')
    parsdata.ntot = [parsdata.ntot; length(pars)]; % Total number of...
        % ...aggregates
else
    parsdata.ntot = [parsdata.ntot; size(pars.n, 1)];
end

if isempty(parsdata.vf)
    if isa(pars, 'AGG')
        pp = AGG.COMPILEPP(pars);
    else
        pp = cell2mat(pars.pp);
    end
    parsdata.vf = pi * sum(pp(:,2) .^ 3) / (6 * params_ud.Value(1) *...
        params_ud.Value(2) * params_ud.Value(3));
end

if isempty(parsdata.beta)
    parsdata.beta = NaN;
else
    kk = length(parsdata.ii);
    parsdata.beta = [parsdata.beta; -2 * (parsdata.ntot(kk) -...
        parsdata.ntot(kk - 1)) /(parsdata.t(kk) - parsdata.t(kk - 1)) /...
        (parsdata.ntot(kk) .^ 2)]; % Collision frequency
end

[dn_dlogdv, dv_discrete] = TRANSP.DISTRIBUTION(ones(length(dv),1), dv);
parsdata.dn_dlogdv = [parsdata.dn_dlogdv; {[dn_dlogdv, dv_discrete]}];
    % Number size distribution

dm_dlogdv = TRANSP.DISTRIBUTION(m, dv);
parsdata.dm_dlogdv = [parsdata.dm_dlogdv, {[dm_dlogdv, dv_discrete]}];
    % Mass size distribution

dg_dpp = dg ./ dpp(:,1); % Aggregate to primary particles size ratio
parsdata.dg_dpp = [parsdata.dg_dpp, {dg_dpp}];
parsdata.npp = [parsdata.npp, {npp}]; % Primaries number distribution
psfit = fit(dg_dpp, npp, 'power1'); % A without-intercept power series...
    % ...fit to the number vs. size ratio variations
parsdata.df = [parsdata.df; psfit.b]; % Fractal dimension
parsdata.kf = [parsdata.kf; psfit.a]; % Fractal prefactor

end
