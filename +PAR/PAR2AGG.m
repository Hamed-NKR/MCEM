function agg = PAR2AGG(par)
% "PAR2AGG" Convert from a PAR structure to an array of AGG objects.
%
%  Original author: Timothy Sipkens, 06-2021
% ----------------------------------------------------------------------- %
% 
% Input:
%     par: Particle information structure
% ----------------------------------------------------------------------- %
% 
% Outputs:
%     agg: Aggregate class objects
% ----------------------------------------------------------------------- %

nagg = length(par.pp);  % number of aggregates

agg = AGG();

for ii=1:nagg
    pp = struct();
    
    pp.r = par.pp{ii}(:, 3:end);
    pp.dp = par.pp{ii}(:, 2);
    pp.id = par.pp{ii}(:, 1);
    
    agg(ii, :) = AGG(pp);
    
    agg(ii).d = par.d(ii, :);
    agg(ii).v = par.v(ii, :);
    
    % TO DO: Mass is not correct.
    agg(ii).m = par.m(ii);                  % ~ mass
    agg(ii).rho = par.rho(ii);              % ~ effective density
    agg(ii).delt = par.delt(ii);        % ~ motion time-step
    agg(ii).tau = par.tau(ii);          % ~ relaxation time
    agg(ii).f = par.f(ii);              % ~ friction factor
    agg(ii).diff = par.diff(ii);        % ~ diffusivity
    agg(ii).lambda = par.lambda(ii);    % ~ diffusive mean free path
    agg(ii).kn = par.kn(ii,:);            % Knudsen number (both kinetic and diffusive)
    agg(ii).nnl = par.nnl{ii};          % ~ nearest neighbor list
end

end

