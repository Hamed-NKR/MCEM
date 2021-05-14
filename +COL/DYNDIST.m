function dist_norm = DYNDIST(r_par, v_par)
% "DYNDIST" computes the minimum center-to-center distance achievable
% between two moving spherical particles.

% Inputs are the particles' location, velocity and size (each having two...
% ...rows corresponding to the pair particles).

v_rel = v_par(2,:) - v_par(1,:);  % relative velocity of two particles
r_rel = r_par(2,:) - r_par(1,:);  % relative distance

dist_norm = norm(cross(r_rel, v_rel)) / norm(v_rel);  % projection distance

end

