function agg = AGG_MOBIL(agg, fl, params_const)
% "AGG_MOBIL" calculates the transport properties of aggregates.
% ----------------------------------------------------------------------- %

% Inputs/output:
    % agg: Aggregate information structure containting its geometrical...
        % ...and diffusive properties
    % fl: Fluid info structure
    % params_const: Problem's table of constant physical properties    
% ----------------------------------------------------------------------- %

n_agg = size(agg.n,1); % Total number of aggregates

agg.kn = zeros(n_agg,2);
agg.kn(:,1) = (2 * fl.lambda) ./ (agg.d); % Kinetic (momentum) Knudsen..
    % ...number
alpha = 1.254; 
beta = 0.4;
gamma = 1.1;
cc = 1 + (agg.kn(:,1)) .* (alpha + beta .* exp(-gamma ./ (agg.kn(:,1))));
    % Cunningham correction factor

rho_bc = params_const.Value(1); % Black Carbon bulk density (kg/m3)
agg.rho = rho_bc; % Effective density
agg.tau = rho_bc .* ((agg.d).^2) .* cc ./ (18 .* fl.mu); % Response...
    % ...(relaxation) time

kb = params_const.Value(3); % Boltzmann's constant (j/k)
agg.f = (3 * pi * fl.mu) .* (agg.d) ./ cc; % Friction factor
agg.diff = (kb * (fl.temp)) ./ (agg.f) ; % Particle diffusivity (m2/s)

for i = 1 : n_agg
    agg.m = (params_const.Value(1) * pi / 6) * sum((agg.pp{i}(:,2)).^3);
        % Aggregate mass (kg)
end

agg.lambda = sqrt((agg.m) .* kb .* (fl.temp)) ./ (agg.f); % Diffusive...
    % ...mean free path

agg.kn(:,2) = (2 * agg.lambda) ./ (agg.d); % Diffusive Knudsen number

% Computing marching timestep
agg.delt = ((agg.f) .* ((agg.d).^2)) ./ (6 .* kb .* fl.temp);

end
