function ovr_pars = OVR(r_pars, d_pars)
% "OVR" checks whether two sets of particle pairs overlap.
% ----------------------------------------------------------------------- %
% 
% Inputs:
%     r_pars: The center location of the particle set pairs (an N*6...
%         ...array as [r_pars1_x, r_pars1_y, r_pars1_z,...
%         ...r_pars2_x, r_pars2_y, r_pars2_z])
%     d_pars: The the particle diameters (N*2: [d_pars1, d_pars2])
% ----------------------------------------------------------------------- %
% 
% Output:
%     ovr_pars: A logical variable indicating whether the particles...
%         ...overlap (1 --> overlapping, 0 --> non-overlapping)
% ----------------------------------------------------------------------- %

if ~ isempty(find(d_pars <= 0,1)) % Check the diameters to be positive
    error("invalid particle diameters!")
end

n_par = size(d_pars,1);
ovr_pars = zeros(n_par,1); % particles initially assumed not overlapping

dist1 = sqrt(sum((r_pars(:,1:3) - r_pars(:,4:6)).^2, 2)); % Particles...
    % ...center-to-center distances
dist2 = (d_pars(:,1) + d_pars(:,2)) / 2; % Overlapping criterion
ovr_pars(dist2 - dist1 > 0) = 1; % particles overlap
 
end

