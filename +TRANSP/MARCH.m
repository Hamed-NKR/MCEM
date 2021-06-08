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

% Compile primary particles across multiple aggregates.
if isa(par, 'AGG')
    n_par = size(par, 1);
else
    n_par = size(par.n, 1); % Total number of particles
end

kb = params_const.Value(3); % Boltzmann's constant (j/k)

% Copy properties locally.
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
% Velocity march
rv_dot_rv = (3 .* kb .* fl.temp ./ m) .* (1 - (var_march.^(-2)));
par_v_new = v .* (var_march.^(-1)) + sqrt(rv_dot_rv ./ 3) .*...
    [randn(n_par, 1), randn(n_par, 1), randn(n_par, 1)];
% position march
rr_dot_rr = (6 .* m .* kb .* fl.temp ./ (f.^2)) .*...
    (f .* delt ./ m - 2 .* ((1 - (var_march.^(-1))) ./...
    (1 + (var_march.^(-1)))));
dr0 = (m ./ f) .* (par_v_new + v) .*...
    ((1 - (var_march.^(-1))) ./ (1 + (var_march.^(-1)))) +...
    sqrt(rr_dot_rr ./ 3) .* [randn(n_par, 1), randn(n_par, 1),...
    randn(n_par, 1)];  % change in position


% Updating the particle structure
if isa(par, 'AGG')
    for ii=1:length(par)  % update velocity
        par(ii).v = par_v_new(ii, :);
    end
    par = par.TRANSLATE(dr0 ./ z_par);
else
    % Updating the primary particle locations
    pp = cell2mat(par.pp);
    pp(:,3:5) = pp(:,3:5) + repelem(dr0, [par.n], 1);
    
    par.v = par_v_new;
    par.r = par.r + dr0 ./ z_par; % Interpolating the...
        % ...displacement for the baseline timestep
    par.pp = mat2cell(pp, par.n);
end

end

