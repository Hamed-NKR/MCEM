function [par, delt_base] = MARCH(par, fl, params_const)
% "MARCH" solves the equation of motion of particles and finds their new...
    % ...location and velocity alomg with their marching time scale.

% Note: For more detail on the equations used here, see a tutorial on...
    % ...on Langevin Dynamics modeling of aerosols by Suresh &...
    % ...Gopalakrishnan (2021).
% ----------------------------------------------------------------------- %

% Input/Outputs:
    % par: Particle information structure
    % fl: Fluid info structure
    % params_const: Problem's table of constant physical properties
    % delt_base: The universal marching times-step for the poulation of...
        % ...particles
% ----------------------------------------------------------------------- %

kb = params_const.Value(3); % Boltzmann's constant (j/k)
n_par = size(par.n,1); % Total number of particles

delt_base = min(par.delt); % The baseline time step chosen as the one...
    % ...for the smallest particle
z_par = ceil(par.delt ./ delt_base); % Integer marching coefficients
par.delt = delt_base .* z_par; % This is to avoid numerical instabilities.

% Solving the equation of motion
var_march = exp(par.f .* par.delt ./ par.m);
% Velocity march
rv_dot_rv = (3 .* kb .* fl.temp ./ par.m) .* (1 - (var_march.^(-2)));
par_v_new = par.v .* (var_march.^(-1)) + sqrt(rv_dot_rv ./ 3) .*...
    [randn(n_par, 1), randn(n_par, 1), randn(n_par, 1)];
% position march
rr_dot_rr = (6 .* par.m .* kb .* fl.temp ./ ((par.f).^2)) .*...
    (par.f .* par.delt ./ par.m - 2 .* ((1 - (var_march.^(-1))) ./...
    (1 + (var_march.^(-1)))));
par_r_new = par.r + (par.m ./ (par.f)) .* (par_v_new + par.v) .*...
    ((1 - (var_march.^(-1))) ./ (1 + (var_march.^(-1)))) +...
    sqrt(rr_dot_rr ./ 3) .* [randn(n_par, 1), randn(n_par, 1),...
    randn(n_par, 1)]; 

% Updating the primary particle locations
pp = cell2mat(par.pp);
pp(:,3:5) = pp(:,3:5) + repelem(par_r_new - par.r, par.n, 1);

% Updating the particle structure
par.v = par_v_new;
par.r = par.r + (par_r_new - par.r) ./ z_par; % Interpolating the...
    % ...displacement for the baseline timestep
par.pp = mat2cell(pp, par.n);

end

