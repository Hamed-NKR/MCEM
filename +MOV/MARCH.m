function [r_par_new, v_par_new] = MARCH(r_par_old, v_par_old,...
    d_par, rho_par, tau_par, temp_f, del_t)
% This function solves for equation of motion of a spherical particle...
% and gives the new location and velocity of the particle as well as...
% its relaxation time.

% Inputs are particle initial position and velociy, size, density, and...
% relaxation time.

% Computing random position and velocity components; for more detail,...
% see Heine & Pratsinis, 2007, "Brownian Coagulation at High Concentration"
k_b = 1.381e-23; % Boltzmann's constant
m_par = rho_par * pi .* (d_par.^3) ./ 6; % Particle mass (kg)
sig_v2 = (k_b * temp_f ./ m_par) .* (1 - exp(-2 .* del_t ./ tau_par));
sig_vr = (k_b * temp_f ./ m_par) .* ((1 - exp(- del_t ./ tau_par)).^2) ./ tau_par;
sig_r2 = (k_b * temp_f ./ m_par) .* ((2 .* del_t ./ tau_par) - 3 ...
    + 4 .* exp(- del_t ./ tau_par) - exp(-2 .* del_t ./ tau_par)) ./ (tau_par.^2);
y = randn(size(d_par,1), 6); % Gaussian-distributed independent random numbers with...
% ...zero mean and variance of unity

% TO DO: There is a bug here, as there is the sqrt() of a negative number.
r_par_rand = (sig_vr./sqrt(sig_v2)) .* y(:,1:3) + ...
    sqrt(sig_r2 - (sig_vr.^2)./sig_v2) .* y(:,4:6); % Random position component
v_par_rand = sqrt(sig_v2) .* y(:,1:3); % Random velocity component

% Finding the new position and velocity vectors
v_par_new = v_par_rand + v_par_old .* exp(-del_t ./ tau_par);
r_par_new = r_par_rand + r_par_old + v_par_old .* tau_par .* ...
    (1 - exp(-del_t ./ tau_par));

end

