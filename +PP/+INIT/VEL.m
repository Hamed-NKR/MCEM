function pp_v = VEL(pp_d,flm_temp)
% This function gives size dependent initial velocities to the...
% primary particles along random orientations.

% Inputs are primary particle size array (mean+std), and flame temperature.

n_pp = size(pp_d,1); % Number of primaries
pp_v = zeros(n_pp,3); % Initialization of the velocity array

rho_bc = 1.8 * 10^3; % Black Carbon's bulk density (kg/m3)
k_b = 1.381 * 10^(-23); % Boltzmann's constant (j/k)
m_pp = (pi *  rho_bc / 6) .* (pp_d.^3); % mass of primaries
% Maxwell-Boltzmann's mean kinetic speed (m/s)
v_mb_mean_pp = ((8 * k_b * flm_temp / pi) ./ m_pp).^(1/2);

% Random spherical angles for initial movements of the primaries  
theta_pp = pi * rand(n_pp,1);
phi_pp = 2 * pi * rand(n_pp,1);

% Assigning random movements to the primaries
pp_v(:,1) = v_mb_mean_pp .* cos(theta_pp) .* sin(phi_pp);
pp_v(:,2) = v_mb_mean_pp .* sin(theta_pp) .* sin(phi_pp);
pp_v(:,3) = v_mb_mean_pp .* cos(phi_pp);

end

