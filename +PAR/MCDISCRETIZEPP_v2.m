function [r, c, n_mc] = MCDISCRETIZEPP_v2(D, CC, n_mc0)
% "MCDISCRETIZEPP" generates a uniformly random (Monte Carlo) set of...
%   ...lattice points within a spherical domain (e.g., primary particle).
% ----------------------------------------------------------------------- %
% 
% Inputs:
%   D: Particle (spherical domain) diameter
%   CC: Central coordinates of the particle
%   n_mc0: Number of expected Monte Carlo grid points within the particle
% ----------------------------------------------------------------------- %
% 
% Outputs:
%   r: spatial location set of lattice points (an n*3 array)
%   c: concentration of lattice points in the domain
%   n_mc: Number of outputted Monte Carlo grid points within the particle
% ----------------------------------------------------------------------- %

% initialize inputs if not defined

if ~exist('D', 'var') || isempty(D)
    D = 1;
end

if ~exist('CC', 'var') || isempty(CC)
    CC = zeros(3,1);
end

if ~exist('n_mc0', 'var') || isempty(n_mc0)
    n_mc0 = 50;
end

% approximate number of grid points in a cube surrounding the sphere
n_mc = ceil(n_mc0 * 6 / pi);

% generate random points in the cube
r = rand(n_mc,3) - 0.5;

% filter the points to get what falls inside the sphere
r = D * r(sqrt(sum(r.^2, 2)) <= 0.5,:) + CC;

n_mc = size(r,1); % update the number of grid points after filtering

c = 6 * n_mc / (pi * D^3); % get the concentration of points within the sphere

