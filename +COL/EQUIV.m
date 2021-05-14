function [r_cm, v_cm] = EQUIV(r_par, v_par, d_par)
% EQUIV computes center of mass and equivalent velocity for an ensemble of
% spherical particles.

% Inputs are particle locations and sizes.

vol_par = pi .* (d_par.^3) ./ 6; % particle volumes
r_cm = sum(vol_par .* r_par) ./ sum(vol_par); % center of mass
v_cm = sum(vol_par .* v_par) ./ sum(vol_par); % equivalent velocity

end

