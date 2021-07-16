function [normdist, normtime] = DYNDIST(r_pps, v_pps)
% "DYNDIST" computes the minimum center-to-center distance achievable...
%     ...between moving spherical particle pairs.
% 
% Note: See the below link for more info:
%     https://math.stackexchange.com/questions/1905533/find-perpendicular-distance-from-point-to-line-in-3d
% ----------------------------------------------------------------------- %
% 
% Inputs:
%     r_pps: Particle locations, an N*6 array as [r_pp1, r_pp2]
%     v_pps: ~ velocities, ~ [v_pp1, v_pp2]
% ----------------------------------------------------------------------- %
% 
% Output:
%     normdist: The closest distance possible between the pairs
%     normtime: The timestep corresponding to the closet distance
% ----------------------------------------------------------------------- %

v_rel = v_pps(:,4:6) - v_pps(:,1:3);  % Relative velocities of particle...
    % ...pairs
r_rel = r_pps(:,4:6) - r_pps(:,1:3);  % Relative displacements ~

% Calculating the relative normal (projection) distance and time to...
    % ...reach it
vnorm_rel = v_rel ./ sqrt(sum(v_rel.^2, 2)); % Normalized relative velocity
d_tang = dot(-r_rel, vnorm_rel, 2); % Tangential distance to the...
    % ...projected point
normtime = abs(d_tang) ./ sqrt(sum(v_rel.^2, 2)); % Projection time
r_prj = r_pps(:,4:6) + d_tang .* vnorm_rel; % Location of projected points
normdist = sqrt(sum((r_prj - r_pps(:,1:3)).^2, 2)); % Projection distance

end
