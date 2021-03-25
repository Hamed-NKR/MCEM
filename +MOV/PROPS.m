function [tau_par, diff_par, mu_f, lambda_f] = PROPS(par_d, fl)
% This function computes the physical propertied of fluid and particles...
% needed for solving the particles equation of motion.

% Inputs are particle size, fluid temperature and pressure.

% Calculating the fluid slip-related properties based on the...
% kinetic theory of gases
R_u = 8.314; % Universal gas constant (j/mol.k)
N_a = 6.022e23; % Avogadro constant (mol^-1)
mm_f = 28.97e-3; % Air molar mass (kg/mol)
d_kin = (0.79 * 3.64 + 0.21 * 3.46) * (e-10); % Air kinetic diameter...
% ...(an average estimate based on values for N2 and O2)
sig_col = pi * (d_kin^2); % Collision area (m2)
mu_f = 2 * sqrt(mm_f*R_u*(fl.temp)) / (3*sqrt(pi)*N_a*sig_col); % Viscosity
lambda_f = mu_f / (0.499*(fl.p)*sqrt(8*mm_f/(pi*R_u*(fl.temp)))); % Mean free path 

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

end

