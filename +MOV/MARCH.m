function [par, delt_base] = MARCH(par, fl)
% "MARCH" solves for equation of motion of a spherical particle and gives
% the new location and velocity of the particle as well as its marching
% time scale.

% For more detail, see the tutorial proposed by Suresh & Gopalakrishnan
% (2021) on Langevin Dynamics modeling of aerosols.

% Inputs are primary particle, aggregate and fluid structures.

n_par = size(par.d,1); % Total number of particles

kn = (2 * fl.lambda) ./ (par.d); % Knudsen number
alpha = 1.254;
beta = 0.4;
gamma = 1.1;
cc = 1 + kn .* (alpha + beta .* exp(-gamma ./ kn));
% Cunningham correction factor

rho_par = 1.8e3; % Primaries density ~ Black carbon's bulk density
par.tau = rho_par .* ((par.d).^2) .* cc ./ (18 .* fl.mu);
% Particle response (relaxation) time

k_b = 1.381e-23; % Boltzmann's constant
f_par = (3 * pi * fl.mu) .* (par.d) .* cc;
par.diff = (k_b * (fl.temp)) ./ f_par; % Particle diffusivity (m2/s)

m_par = rho_par .* pi .* ((par.d).^3) ./ 6; % Particle mass (kg)
par.lambda = sqrt(m_par.*k_b.*(fl.temp)) ./ f_par; % Particle diffusive...
% ...mean free path

% Computing marching timestep
par.delt = (f_par .* ((par.d).^2)) ./ (6 .* k_b .* fl.temp);
% Preliminary time steps
delt_base = min(par.delt); % baseline time step from the smallest particle
z_par = ceil((par.delt)./delt_base); % Integer marching coefficients
par.delt = delt_base .* z_par; % This is to avoid numerical instabilities.

% Solving the equation of motion
var_march = exp(f_par .* (par.delt) ./m_par);
% Velocity march
rv_dot_rv = (3 .* k_b .* fl.temp ./ m_par) .* (1 - (var_march.^(-2)));
par_v_new = par.v .* (var_march.^(-1)) + sqrt(rv_dot_rv./3) .* ...
    [randn(n_par, 1), randn(n_par, 1), randn(n_par, 1)];
% position march
rr_dot_rr = (6 .* m_par .* k_b .* fl.temp ./ (f_par.^2)) .* ...
    (f_par .* (par.delt) ./m_par - 2 .* ((1 - (var_march.^(-1))) ./ ...
    (1 + (var_march.^(-1)))));
par_r_new = par.r + (m_par ./ f_par) .* (par_v_new + par.v) .* ...
    ((1 - (var_march.^(-1))) ./ (1 + (var_march.^(-1)))) + ...
    sqrt(rr_dot_rr./3) .* [randn(n_par, 1), randn(n_par, 1), ...
    randn(n_par, 1)]; 

% Updating the particle structure
par.v = par_v_new;
par.r = par.r + (par_r_new - par.r) ./ z_par;
% Interpolating the displacement for the baseline timestep

end

