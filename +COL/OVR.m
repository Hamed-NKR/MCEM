function ovr_par = OVR(r_par,d_par)
% "OVR" checks whether two particles overlap.

% r_par is 2*3 array for the center location of two particles.
% d_par is 2*1 vector for the particle diameters.

dist1 = sqrt(sum((r_par(2,:) - r_par(1,:)).^2)); % Particle centers'...
% ...distance
dist2 = (d_par(1) + d_par(2)) / 2; % Overlapping criterion

if ~ isempty(find(d_par <= 0,1)) % Check the diameters to be positive
    error("invalid particle diameters!")
elseif dist1 <= dist2
    ovr_par = 1; % particles overlap
else
    ovr_par = 0; % particles do not overlap
end

