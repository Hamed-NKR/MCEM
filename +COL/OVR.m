function par_ovr = OVR(par_loc,par_diam)
% This function checks whether two particles overlap.

% pp_loc is 2*3 array for the center location of two particles.
% pp_diam is 2*1 vector for the particles' diameters.

dist1 = sqrt(sum((par_loc(1,:) - par_loc(2,:)).^2)); % Particle centers' distance
dist2 = (par_diam(1) + par_diam(2)) / 2; % Overlapping criterion

if ~ isempty(find(par_diam <= 0,1)) % Check the diameters to be positive
    error("invalid particle diameters!")
elseif dist1 <= dist2
    par_ovr = 1; % particles overlap
else
    par_ovr = 0; % particles do not overlap
end

