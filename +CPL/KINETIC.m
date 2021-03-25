function [mu_f, lambda_f] = KINETIC(fl_temp,fl_p)
% This function Calculates the fluid slip-related properties based on...
% ...the kinetic theory of gases

% Inputs are the fluid temperature and pressure.

R_u = 8.314; % Universal gas constant (j/mol.k)
N_a = 6.022e23; % Avogadro constant (mol^-1)
mm_f = 28.97e-3; % Air molar mass (kg/mol)
d_kin = (0.79 * 3.64 + 0.21 * 3.46) * (e-10); % Air kinetic diameter...
% ...(an average estimate based on values for N2 and O2)
sig_col = pi * (d_kin^2); % Collision area (m2)
mu_f = 2 * sqrt(mm_f*R_u*(fl_temp)) / (3*sqrt(pi)*N_a*sig_col); % Viscosity
lambda_f = mu_f / (0.499*(fl_p)*sqrt(8*mm_f/(pi*R_u*(fl_temp)))); % Mean free path 

end

