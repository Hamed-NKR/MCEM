function vec = ROTATE_SIMPLE(vec, r_c, dtheta)
% "ROTATE_SIMPLE" rotates a vector around an arbitrary center.
% ----------------------------------------------------------------------- %
% 
% Inputs/Outputs:
%     vec: The vector considered to rotate
%     r_c: Center of rotation
%     dtheta: Rotation vector (1/ yaw, 2/ pitch, 3/ roll), an N*3 vector
% ----------------------------------------------------------------------- %

n_vec = size(vec, 1); % Number of vectors to be rotated

dtheta = reshape(dtheta, 1, 3, n_vec);
vec = reshape(vec, 1, 3, n_vec);
r_c = reshape(r_c, 1, 3, n_vec);

yaw = [cos(dtheta(1,1,:)), -sin(dtheta(1,1,:)), zeros(1,1,n_vec);...
    sin(dtheta(1,1,:)), cos(dtheta(1,1,:)), zeros(1,1,n_vec);...
    zeros(1,1,n_vec), zeros(1,1,n_vec), ones(1,1,n_vec)];
    % transformation matrix for yaw rotation
pitch = [cos(dtheta(1,2,:)), zeros(1,1,n_vec), sin(dtheta(1,2,:));...
     zeros(1,1,n_vec), ones(1,1,n_vec), zeros(1,1,n_vec);...
     -sin(dtheta(1,2,:)), zeros(1,1,n_vec), cos(dtheta(1,2,:))];
    % transformation matrix for pitch rotation
roll = [ones(1,1,n_vec), zeros(1,1,n_vec), zeros(1,1,n_vec);...
    zeros(1,1,n_vec), cos(dtheta(1,3,:)), -sin(dtheta(1,3,:));...
    zeros(1,1,n_vec), sin(dtheta(1,3,:)), cos(dtheta(1,3,:))];
    % transformation matrix for roll rotation

% The net rotation matrix
rot = pagemtimes(yaw, pitch);
rot = pagemtimes(rot, roll);

% Rotating the vector
vec = permute(pagemtimes(rot,'none', (vec - r_c), 'transpose'), [2 1 3])...
    + r_c;
vec = reshape(vec, n_vec, 3, 1);

end

