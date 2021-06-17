function par = MOBIL(par, fl, params_const)
% "MOBIL" calculates the transport properties of aggregates.
% ----------------------------------------------------------------------- %

% Inputs/output:
    % par: Aggregates information structure containting their...
        % ...geometrical and diffusive properties
    % fl: Fluid info structure
    % params_const: Problem's table of constant physical properties    
% ----------------------------------------------------------------------- %

n_par = size(par.n,1); % Total number of aggregates

par.kn = zeros(n_par,2);
par.kn(:,1) = (2 * fl.lambda) ./ (par.d); % Kinetic (momentum) Knudsen..
    % ...number
alpha = 1.254; 
beta = 0.4;
gamma = 1.1;
cc = 1 + (par.kn(:,1)) .* (alpha + beta .* exp(-gamma ./ (par.kn(:,1))));
    % Cunningham correction factor

rho_bc = params_const.Value(1); % Black Carbon bulk density (kg/m3)
par.rho = rho_bc .* ones(n_par,1); % Effective density
par.tau = rho_bc .* ((par.d).^2) .* cc ./ (18 .* fl.mu); % Response...
    % ...(relaxation) time

kb = params_const.Value(3); % Boltzmann's constant (j/k)
par.f = (3 * pi * fl.mu) .* (par.d) ./ cc; % Friction factor
par.diff = (kb * (fl.temp)) ./ (par.f) ; % Particle diffusivity (m2/s)

if isempty(par.m)
    par.m = ones(n_par,1); % Initializing aggregate mass array
    for i = 1 : n_par
        par.m(i) = (params_const.Value(1) * pi / 6) *...
            sum((par.pp{i}(:,2)).^3); % Aggregate mass (kg)
    end
end

par.lambda = sqrt((par.m) .* kb .* (fl.temp)) ./ (par.f); % Diffusive...
    % ...mean free path

par.kn(:,2) = (2 * par.lambda) ./ (par.d); % Diffusive Knudsen number

% Computing marching timestep
par.delt = ((par.f) .* ((par.d).^2)) ./ (6 .* kb .* fl.temp);

end
