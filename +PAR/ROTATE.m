function par_pp = ROTATE(par_pp, par_r, par_n, dtheta)
% "ROTATE" rotates the primary particle locations of an aggregate...
%     ...with respect to their center of mass.
% ----------------------------------------------------------------------- %
% 
% Inputs/Outputs:
%     par_pp: Primary particles info cell array
%     par_r: Center of mass of (independent) particles
%     par_n: Number of primaries within each aggregate
%     dtheta: Rotation vector (1/ yaw, 2/ pitch, 3/ roll), an N*3 vector
% ----------------------------------------------------------------------- %

n_par = size(par_n, 1); % Total number of aggregates

dtheta = reshape(dtheta, 1, 3, n_par);

yaw = [cos(dtheta(1,1,:)), -sin(dtheta(1,1,:)), zeros(1,1,n_par);...
    sin(dtheta(1,1,:)), cos(dtheta(1,1,:)), zeros(1,1,n_par);...
    zeros(1,1,n_par), zeros(1,1,n_par), ones(1,1,n_par)];
    % transformation matrix for yaw rotation
pitch = [cos(dtheta(1,2,:)), zeros(1,1,n_par), sin(dtheta(1,2,:));...
     zeros(1,1,n_par), ones(1,1,n_par), zeros(1,1,n_par);...
     -sin(dtheta(1,2,:)), zeros(1,1,n_par), cos(dtheta(1,2,:))];
    % transformation matrix for pitch rotation
roll = [ones(1,1,n_par), zeros(1,1,n_par), zeros(1,1,n_par);...
    zeros(1,1,n_par), cos(dtheta(1,3,:)), -sin(dtheta(1,3,:));...
    zeros(1,1,n_par), sin(dtheta(1,3,:)), cos(dtheta(1,3,:))];
    % transformation matrix for roll rotation

% The net rotation matrix
rot = pagemtimes(yaw, pitch);
rot = pagemtimes(rot, roll);

% Obtaining the post-rotation coordinates of primaries
for i = 1 : n_par
    par_pp{i}(:,3:5) = (rot(:,:,i) * (par_pp{i}(:,3:5) -...
        repmat(par_r(i,:), par_n(i), 1))')' +...
        repmat(par_r(i,:), par_n(i), 1);
end

end

