function aggs = PAR2AGG(pars)
% "PAR2AGG" Convert from a PAR structure to an array of AGG objects.
%
% Original author: Timothy Sipkens, 06-2021
% Revised by: Hamed Nikookar, 07-2021
% ----------------------------------------------------------------------- %
% 
% Input:
%     pars: A global structure containing particles information 
% ----------------------------------------------------------------------- %
% 
% Output:
%     aggs: Objects of aggregate class
% ----------------------------------------------------------------------- %

n_agg = size(pars.n,1);  % Number of aggregates

aggs = AGG(); % Initializng the aggregate class

for i = 1 : n_agg
    
    % Extracting the primary particles data
    pp = struct();
    pp.id = pars.pp{i}(:,1);
    pp.d = pars.pp{i}(:,2);
    pp.r = pars.pp{i}(:,3:5);
    
    % Generating aggregate objects
    aggs(i) = AGG(pp);
    
    % Assigning the aggregates physical properties
    aggs(i).dv = pars.dv(i);               % Volumetric diameter
    aggs(i).dg = pars.dg(i);               % Gyration diameter
%     aggs(i).dm = pars.dm(i);               % Mobility diameter
%     aggs(i).da = pars.da(i);               % Aerodynamic diameter
    aggs(i).dpp = pars.dpp(i,:);           % Maximum extent diameter
    aggs(i).dmax = pars.dmax(i);           % Maximum extent diameter
    
    aggs(i).r = pars.r(i,:);               % Location
    aggs(i).v = pars.v(i,:);               % Velocity
    
    aggs(i).m = pars.m(i);                 % Mass
    aggs(i).rho = pars.rho(i);             % Effective density
    aggs(i).delt = pars.delt(i);           % Motion time-step
    aggs(i).tau = pars.tau(i);             % Relaxation time
    aggs(i).f = pars.f(i);                 % Friction factor
    aggs(i).diff = pars.diff(i);           % Diffusivity
    aggs(i).lambda = pars.lambda(i);       % Diffusive mean free path
    aggs(i).kn_kin = pars.kn_kin(i,:);     % Kinetic Knudsen number
    aggs(i).kn_diff = pars.kn_diff(i,:);   % Diffusive Knudsen number
    aggs(i).nnl = pars.nnl{i};             % Nearest neighbor list
    
end

end

