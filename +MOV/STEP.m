function [delt_par, tau_par] = STEP(par_d, fl)
% This function computes the marching time step of particles...
% needed for solving their equation of motion.

% Inputs are particle size and fluid structure.


% Calculating the particle slip-related properties
kn = 2 .* lambda_f ./ par_d; % Knudsen number
alpha = 1.254;
beta = 0.4;
gamma = 1.1;
cc = 1 + kn .* (alpha + beta .* exp(-gamma ./ kn)); % Cunningham correction factor

rho_par = 1.8e3; % Primaries density ~ Black carbon's bulk density
tau_par = rho_par .* (par_d.^2) .* cc ./ (18 .* mu_f);
% Particle response (relaxation) time

k_b = 1.381e-23; % Boltzmann's constant
f_par = (3 * pi * mu_f) .* par_d;
diff_par = (k_b * (fl.temp)) ./ f_par; % Particle diffusivity (m2/s)
delt_par = ((0.02 .* par_d).^2) ./ (2 .* diff_par); % Time step

end

