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
    parsdata = struct('ii', [], 't', [], 'n_tot', [], 'vf', [],...
        'beta', [], 'dn_dlogdv', [], 'dm_dlogdv', []);
end

parsdata.ii = [parsdata.ii; k]; % Adding iteration index
parsdata.t = [parsdata.t; time(k)]; % Time

if isa(pars, 'AGG')
    parsdata.n_tot = [parsdata.n_tot; length(pars)]; % Total number of...
        % ...aggregates
else
    parsdata.n_tot = [parsdata.n_tot; size(pars.n, 1)];
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
    parsdata.beta = [parsdata.beta; -2 * (parsdata.n_tot(kk) -...
        parsdata.n_tot(kk - 1)) /(parsdata.t(kk) - parsdata.t(kk - 1)) /...
        (parsdata.n_tot(kk) .^ 2)]; % Collision frequency
end

dv = cat(1, pars.dv);

[dn_dlogdv, dv_discrete] = TRANSP.DISTRIBUTION(ones(length(dv),1), dv);
parsdata.dn_dlogdv = [parsdata.dn_dlogdv; {[dn_dlogdv, dv_discrete]}];
    % Number size distribution

m = cat(1, pars.m);
dm_dlogdv = TRANSP.DISTRIBUTION(m, dv);
parsdata.dm_dlogdv = [parsdata.dm_dlogdv, {[dm_dlogdv, dv_discrete]}];
    % Mass size distribution

end
