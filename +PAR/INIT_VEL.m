function par_v = INIT_VEL(pp, n_pp ,fl_temp, params_const)
% "INIT_VEL" gives size dependent initial velocities to the primary...
%     ...particles along random orientations.
% ----------------------------------------------------------------------- %
% 
% Inputs:
%     par_pp: Primary particle infomration cell array
%     par_n: Number distribution of primaries
%     fl_temp: Flame temperature
%     params_const: Problem's table of constant physical properties
% ----------------------------------------------------------------------- %
% 
% Output:
%     par_v: Primary particle velocities
% ----------------------------------------------------------------------- %

rho_bc = params_const.Value(1); % Black Carbon's bulk density (kg/m3)
kb = params_const.Value(3); % Boltzmann's constant (j/k)

n_par = size(n_pp,1); % Total number of (independent) particles
pp = cell2mat(pp);
m_pp = (pi *  rho_bc / 6) .* ((pp(:,2)).^3); % mass of primaries
m_pp = mat2cell(m_pp, n_pp);
m_par = zeros(n_par,1);
for i = 1 : n_par
    m_par(i) = sum(m_pp{i}); % Total mass of (independent) particles
end 

vmb_mean_par = ((8 * kb * fl_temp / pi) ./ m_par).^(1/2); % Maxwell-...
    % ...Boltzmann's mean velocity of the (independent) particles

% Random spherical angles for initial movements of the primaries  
theta_pp = pi * rand(n_par,1);
phi_pp = 2 * pi * rand(n_par,1);

par_v = zeros(n_par,3); % Initialization of the velocity array
% Assigning random movements to the primaries
par_v(:,1) = vmb_mean_par .* sin(theta_pp) .* cos(phi_pp);
par_v(:,2) = vmb_mean_par .* sin(theta_pp) .* sin(phi_pp);
par_v(:,3) = vmb_mean_par .* cos(theta_pp);

end

