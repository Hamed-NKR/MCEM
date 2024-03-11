function r = CUBEDISCRETIZEPP(D, C, n_rsl)
% "CUBEDISCRETIZEPP" generates a 3d structured set of lattice points...
%   ...within a spherical domain.
% ----------------------------------------------------------------------- %
% 
% Inputs:
%   D: Particle (sphere) diameter
%   C: Central coordinates of the lattice set
%   n_rsl: resolution variable (number of increments in each dimension)
% ----------------------------------------------------------------------- %
% 
% Outputs:
%   r: spatial location of lattice points (an n*3 array)
% ----------------------------------------------------------------------- %

% initialize inputs if not defined

if ~exist('D', 'var') || isempty(D)
    D = 1;
end

if ~exist('C', 'var') || isempty(C)
    C = zeros(3,1);
end

if ~exist('n_rsl', 'var') || isempty(n_rsl)
    n_rsl = 10;
end

% initialize coordinates
rx = (0 : D / n_rsl : D) - D / 2;
ry = (0 : D / n_rsl : D) - D / 2;
rz = (0 : D / n_rsl : D) - D / 2;
r = zeros((n_rsl + 1)^3, 3);

% make a cubic structured grid
for i = 1 : (n_rsl + 1)
    for j = 1 : (n_rsl + 1)
        for k = 1 : (n_rsl + 1)
            r(i + (j-1) * (n_rsl + 1) + (k-1) * (n_rsl +1)^2, :) =...
                [rx(i), ry(j), rz(k)];
        end
    end
end

C = C(1:3);
C = reshape(C, [1,3]);
% filter the grid (to points within the particle) and move grid to the particle center
r = r(sqrt(sum(r.^2, 2)) <= D/2,:) + C;





