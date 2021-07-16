function pars = MOBIL(pars, fl, params_const)
% "MOBIL" calculates the transport properties of aggregates.
% ----------------------------------------------------------------------- %
% 
% Inputs/output:
%     pars: Aggregates information structure containting their...
%         ...geometrical and diffusive properties
%     fl: Fluid info structure
%     params_const: Problem's table of constant physical properties    
% ----------------------------------------------------------------------- %

% Total number of (independent) particles
if isa(pars, 'AGG')
    n_par = length(pars);
else
    n_par = size(pars.n, 1);
end

% Compiling/copying properties locally
dv = cat(1, pars.dv);

kn_kin = (2 * fl.lambda) ./ dv; % Kinetic (momentum) Knudsen number (-)
alpha = 1.254; 
beta = 0.4;
gamma = 1.1;
cc = 1 + (kn_kin(:,1)) .* (alpha + beta .* exp(-gamma ./ (kn_kin(:,1))));
    % Cunningham correction factor (-)

rho_bc = params_const.Value(1); % Black Carbon bulk density (kg/m3)
rho = rho_bc .* ones(n_par,1); % Effective density (kg/m3)
tau = rho .* (dv.^2) .* cc ./ (18 * fl.mu); % Response (relaxation) time (s)
kb = params_const.Value(3); % Boltzmann's constant (j/k)
f = (3 * pi * fl.mu) .* (dv) ./ cc; % Friction factor (-)
delt = f .* (dv.^2) ./ (6 * kb * fl.temp); % Marching timestep
diff = (kb * (fl.temp)) ./ f; % Particle diffusivity (m2/s)

% Aggregate mass (kg)
m = ones(n_par,1);
for i = 1 : n_par
    if isa(pars, 'AGG')
        m(i) = (params_const.Value(1) * pi / 6) * sum((pars(i).pp.d(:)).^3);
    else
        m(i) = (params_const.Value(1) * pi / 6) * sum((pars.pp{i}(:,2)).^3);
    end
end

lambda = sqrt(m .* kb .* (fl.temp)) ./ f; % Diffusive mean free path (m)
kn_diff = 2 * lambda ./ dv; % Diffusive Knudsen number (-)

% Updating the transport properties for the structure/class of aggregates
if isa(pars, 'AGG')
    for i = 1 : n_par
        pars(i).rho = rho(i);
        pars(i).m = m(i);
        pars(i).tau = tau(i);
        pars(i).f = f(i);
        pars(i).delt = delt(i);
        pars(i).diff = diff(i);
        pars(i).lambda = lambda(i);
        pars(i).kn_kin = kn_kin(i);
        pars(i).kn_diff = kn_diff(i);
    end

else
    pars.rho = rho;
    pars.m = m;
    pars.tau = tau;
    pars.f = f;
    pars.delt = delt;
    pars.diff = diff;
    pars.lambda = lambda;
    pars.kn_kin = kn_kin;
    pars.kn_diff = kn_diff;
end

end
