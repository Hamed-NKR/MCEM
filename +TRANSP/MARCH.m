function [par, delt_base] = MARCH(par, fl, params_const)
% "MARCH" solves the equation of motion of particles and finds their new...
%   ...location and velocity alomg with their marching time scale.
%   
% Note: For more detail on the equations used here, see a tutorial on...
%   ...on Langevin Dynamics modeling of aerosols by Suresh &...
%   ...Gopalakrishnan (2021).
% ----------------------------------------------------------------------- %
% 
% Input/Outputs:
%   par: Particle information structure
%   fl: Fluid info structure
%   params_const: Problem's table of constant physical properties
%   delt_base: The universal marching times-step for the poulation of...
%   ...particles
% ----------------------------------------------------------------------- %

% Total number of (independent) particles
if isa(par, 'AGG')
    n_par = size(par, 1);
else
    n_par = size(par.n, 1);
end

kb = params_const.Value(3); % Boltzmann's constant (j/k)

% Compiling/copying properties locally
m = cat(1, par.m);
v = cat(1, par.v);
f = cat(1, par.f);
delt0 = cat(1, par.delt);

delt_base = min(delt0); % The baseline time step chosen as the one...
    % ...for the smallest particle
z_par = ceil(delt0 ./ delt_base); % Integer marching coefficients
delt = delt_base .* z_par; % This is to avoid numerical instabilities.

% Solving the equation of motion
var_march = exp(f .* delt ./ m);
% Finding the velocity march
rv_dot_rv = (3 .* kb .* fl.temp ./ m) .* (1 - (var_march.^(-2)));
par_v_new = v .* (var_march.^(-1)) + sqrt(rv_dot_rv ./ 3) .*...
    [randn(n_par, 1), randn(n_par, 1), randn(n_par, 1)];
% Finding the position march
rr_dot_rr = (6 .* m .* kb .* fl.temp ./ (f.^2)) .*...
    (f .* delt ./ m - 2 .* ((1 - (var_march.^(-1))) ./...
    (1 + (var_march.^(-1)))));
dr = (m ./ f) .* (par_v_new + v) .*...
    ((1 - (var_march.^(-1))) ./ (1 + (var_march.^(-1)))) +...
    sqrt(rr_dot_rr ./ 3) .* [randn(n_par, 1), randn(n_par, 1),...
    randn(n_par, 1)];  % Translation vectors

if isa(par, 'AGG') % Updating aggregate class objects
    for i = 1 : n_par
        par(i).v = par_v_new(i, :);
    end
    par = par.TRANSLATE(dr ./ z_par);

else % Updating the particle structure
    % ~ the primary particle and aggregate locations
    [par.pp, par.r] = PAR.TRANSLATE(par.pp, par.r, par.n, dr ./ z_par);
        % The displacements are interpoldated for the baseline timestep
    par.v = par_v_new; % ~ the velocities

end

end

