function [r_par_new, v_par_new, delt_base] = MARCH(par_old, fl)
% This function solves for equation of motion of a spherical particle...
% and gives the new location and velocity of the particle as well as...
% its marching time scale.

% For more detail, see Heine & Pratsinis, 2007,...
% "Brownian Coagulation at High Concentration"

% Inputs are particle and fluid structures.

rho_par = 1.8e3; % Primaries density ~ Black carbon's bulk density

% Primaries characteristic time and dissusivity
[tau_par, diff_par] = MOV.SLIP(par_old.d,fl); 

% Computing marching timestep
delt_par = ((1e-4 .* (par_old.d)).^2) ./ (2 .* diff_par);
% Preliminary time steps
delt_base = min(delt_par); % baseline time step from the smallest particle
z_par = ceil(delt_par./delt_base); % Integer marching coefficients
delt_par = delt_base .* z_par; % This is to avoid numerical instabilities.

% Computing random position and velocity components.
k_b = 1.381e-23; % Boltzmann's constant
m_par = rho_par * pi .* ((par_old.d).^3) ./ 6; % Particle mass (kg)
sig_v2 = (k_b * (fl.temp) ./ m_par) .* (1 - exp(-2 .* delt_par ./ tau_par));
sig_vr = (k_b * (fl.temp) ./ m_par) .* ((1 - exp(- delt_par ./ tau_par)).^2)...
    ./ tau_par;
sig_r2 = (k_b * (fl.temp) ./ m_par) .* ((2 .* delt_par ./ tau_par) - 3 ...
    + 4 .* exp(- delt_par ./ tau_par) - exp(-2 .* delt_par ./ tau_par)) ...
    ./ (tau_par.^2);
y = randn(size(par_old.d,1), 6); % Gaussian-distributed independent...
% ...random numbers with zero mean and variance of unity

% TO DO: There is a bug here, as there is the sqrt() of a negative number.
r_par_rand = (sig_vr./sqrt(sig_v2)) .* y(:,1:3) + ...
    sqrt(sig_r2 - (sig_vr.^2)./sig_v2) .* y(:,4:6); % Random position component
v_par_rand = sqrt(sig_v2) .* y(:,1:3); % Random velocity component

% Finding the new position and velocity vectors
v_par_new = (v_par_rand + (par_old.v) .* exp(-delt_par ./ tau_par));
r_par_new = par_old.r + (r_par_rand + (par_old.v) .* tau_par .* ...
    (1 - exp(-delt_par ./ tau_par))) ./ z_par;
v_par_new = par_old.v + (v_par_new - (par_old.v)) ./ z_par; 

end

