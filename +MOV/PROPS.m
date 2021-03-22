function [tau_par, mu_f, lambda_f] = PROPS(d_par, rho_par, temp_f, p_f)
% This function computes the physical propertied of fluid and particles...
% needed for solving the particles equation of motion.

% Inputs are particle size and density, and fluid temperature and pressure.

% Calculating the fluid slip-related properties based on the...
% kinetic theory of gases
R_u = 8.314; % Universal gas constant (j/mol.k)
N_a = 6.022 * 10^(23); % Avogadro constant (mol^-1)
mm_f = 28.97 * 10^(-3); % Air molar mass (kg/mol)
d_kin = 0.79 * 3.64 * 10^(-10) + 0.21 * 3.46 * 10^(-10); % Air kinetic diameter...
% ...(an average estimate based on values for N2 and O2)
sig_col = pi * (d_kin^2); % Collision area (m2)
mu_f = 2 * sqrt(mm_f*R_u*temp_f) / (3*sqrt(pi)*N_a*sig_col); % Viscosity
lambda_f = mu_f / (0.499*p_f*sqrt(8*mm_f/(pi*R_u*temp_f))); % Mean free path 

% Calculating the particle slip-related properties
kn = 2 * lambda_f / d_par; % Knudsen number
alpha = 1.254;
beta = 0.4;
gamma = 1.1;
cc = 1 + kn * (alpha + beta * exp(-gamma/kn)); % Cunningham correction factor

% Particle response (relaxation) time
tau_par = rho_par * (d_par^2) * cc / (18 * mu_f);

end

