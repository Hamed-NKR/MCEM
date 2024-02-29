function [r, c] = MCDISCRETIZEPP(D, CC, n_mc)
% "MCDISCRETIZEPP" generates a uniformly random (Monte Carlo) set of...
%   ...lattice points within a spherical domain (e.g., primary particle).
% ----------------------------------------------------------------------- %
% 
% Inputs:
%   D: Particle (spherical domain) diameter
%   CC: Central coordinates of the particle
%   n_mc: Number of Monte Carlo grid points within the particle
% ----------------------------------------------------------------------- %
% 
% Outputs:
%   r: spatial location set of lattice points (an n*3 array)
%   c: concentration of lattice points in the domain
% ----------------------------------------------------------------------- %

% initialize inputs if not defined

if ~exist('D', 'var') || isempty(D)
    D = 1;
end

if ~exist('CC', 'var') || isempty(CC)
    CC = zeros(3,1);
end

if ~exist('n_rsl', 'var') || isempty(n_mc)
    n_mc = 50;
end

% generate uniformly random spherical coordinates
r0 = rand(n_mc,3); % columns correspond to rho, theta and phi values respectively
r = zeros(n_mc,3); % initialize coordinates
r(:,1) = (D / 2) * r0(:,1) .* cos(2 * pi * r0(:,2)) .* sin(pi * r0(:,3)) + CC(1); % x = rho * cos(theta) * sin(phi) 
r(:,2) = (D / 2) * r0(:,1) .* sin(2 * pi * r0(:,2)) .* sin(pi * r0(:,3)) + CC(2); % y = rho * sin(theta) * sin(phi) 
r(:,3) = (D / 2) * r0(:,1) .* cos(pi * r0(:,3)) + CC(3); % z = rho * cos(phi)

c = 6 * n_mc / (pi * D^3);

