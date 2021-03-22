function pp_ovr = OVR(pp_loc,pp_diam)
% This function checks whether two primary particles overlap.

% pp_loc is 2*3 array for the center location of two particles.
% pp_diam is 2*1 vector for the particles' diameters.

dist1 = sqrt(sum((pp_loc(1,:) - pp_loc(2,:)).^2)); % Particle centers' distance
dist2 = (pp_diam(1) + pp_diam(2)) / 2; % Overlapping criterion

if ~ isempty(find(pp_diam <= 0,1)) % Check the diameters to be positive
    error("invalid particle diameters!")
elseif dist1 <= dist2
    pp_ovr = 1; % particles overlap
else
    pp_ovr = 0; % particles do not overlap
end

