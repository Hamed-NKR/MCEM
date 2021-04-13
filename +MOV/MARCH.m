function [pp, agg, delt_base] = MARCH(pp, agg, fl)
% This function solves for equation of motion of a spherical particle...
% and gives the new location and velocity of the particle as well as...
% its marching time scale.

% For more detail, see: Suresh, V., & Gopalakrishnan, R. (2021). ...
% Tutorial: Langevin Dynamics methods for aerosol particle trajectory ...
% simulations and collision rate constant modeling. ...
% Journal of Aerosol Science, 155, 105746.

% Inputs are primary particle, aggregate and fluid structures.

n_pp = size(pp.d,1); % Number of primaries
n_agg = size(agg.d,1); % Number of aggregates
n_par = n_pp + n_agg; % Total number of particles
par_d = [pp.d; agg.d];

kn = (2 * fl.lambda) ./ (par_d); % Knudsen number
alpha = 1.254;
beta = 0.4;
gamma = 1.1;
cc = 1 + kn .* (alpha + beta .* exp(-gamma ./ kn)); % Cunningham correction factor

rho_par = 1.8e3; % Primaries density ~ Black carbon's bulk density
par_tau = rho_par .* ((par_d).^2) .* cc ./ (18 .* fl.mu);
% Particle response (relaxation) time

k_b = 1.381e-23; % Boltzmann's constant
f_par = (3 * pi * fl.mu) .* (par_d) .* cc;
par_diff = (k_b * (fl.temp)) ./ f_par; % Particle diffusivity (m2/s)

m_par = rho_par .* pi .* ((par_d).^3) ./ 6; % Particle mass (kg)
par_lambda = sqrt(m_par.*k_b.*(fl.temp)) ./ f_par; % Particle diffusive...
% ...mean free path

% Computing marching timestep
par_delt = (f_par .* ((par_d).^2)) ./ (6 .* k_b .* fl.temp);
% Preliminary time steps
delt_base = min(par_delt); % baseline time step from the smallest particle
z_par = ceil((par_delt)./delt_base); % Integer marching coefficients
par_delt = delt_base .* z_par; % This is to avoid numerical instabilities.

% Solving equation of motion
var_march = exp(f_par .* (par_delt) ./m_par);
par_v_old = [pp.v; agg.v];
par_r_old = [pp.r; agg.r];
% Velocity march
rv_dot_rv = (3 .* k_b .* fl.temp ./ m_par) .* (1 - (var_march.^(-2)));
par_v_new = par_v_old .* (var_march.^(-1)) + sqrt(rv_dot_rv./3) .* ...
    [randn(n_par, 1), randn(n_par, 1), randn(n_par, 1)];
% position march
rr_dot_rr = (6 .* m_par .* k_b .* fl.temp ./ (f_par.^2)) .* ...
    (f_par .* (par_delt) ./m_par - 2 .* ((1 - (var_march.^(-1))) ./ ...
    (1 + (var_march.^(-1)))));
par_r_new = par_r_old + (m_par ./ f_par) .* (par_v_new + par_v_old) .* ...
    ((1 - (var_march.^(-1))) ./ (1 + (var_march.^(-1)))) + ...
    sqrt(rr_dot_rr./3) .* [randn(n_par, 1), randn(n_par, 1), ...
    randn(n_par, 1)]; 

% Interpolating displacement for the baseline timestep and updating the...
% ...particle structure
par_r_new = par_r_old + (par_r_new - par_r_old) ./ z_par;

% Updating the primary particle and aggregate structures
pp.r = par_r_new(1 : n_pp, :);
agg.r = par_r_new(n_pp + 1 : n_pp + n_agg, :);
pp.v = par_v_new(1 : n_pp, :);
agg.v = par_v_new(n_pp + 1 : n_pp + n_agg, :);
pp.tau = par_tau(1 : n_pp);
agg.tau = par_tau(n_pp + 1 : n_pp + n_agg);
pp.diff = par_diff(1 : n_pp);
agg.diff = par_diff(n_pp + 1 : n_pp + n_agg);
pp.lambda = par_lambda(1 : n_pp);
agg.lambda = par_lambda(n_pp + 1 : n_pp + n_agg);
pp.delt = par_delt(1 : n_pp);
agg.delt = par_delt(n_pp + 1 : n_pp + n_agg);

end

