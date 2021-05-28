function r_par_new = PBC(dom_size, r_par_old)
% "PBC" applies periodic boundary condition to the particle movements.

% Inputs are the computational domain size, and particle locations.

n_par = size(r_par_old,1); % Total number of particles
dom_size = repmat(dom_size',n_par,1);
stat_par = r_par_old ./ dom_size; % This status array shows...
% whether the particles are inside the domain.

% Updating position array in the case of domain escape
r_par_new = r_par_old;
if ~ isempty(find(stat_par(:,1) > 1,1)) % checking particle exits from +x
    inds_plus_x = find(stat_par(:,1) > 1); % finding the out of domain...
    % ... particles
    r_par_new(inds_plus_x,1) = mod(r_par_old(inds_plus_x,1), ...
        dom_size(inds_plus_x,1)); % Replacing the outgoing particles with...
    % ...new ones coming from the opposite side.
elseif ~ isempty(find(stat_par(:,1) < 0,1)) % checking particle exits from -x
    inds_minus_x = find(stat_par(:,1) < 0);
    r_par_new(inds_minus_x,1) = dom_size(inds_minus_x,1) - ...
        mod(abs(r_par_old(inds_minus_x,1)), dom_size(inds_minus_x,1));
end

if ~ isempty(find(stat_par(:,2) > 1,1)) % checking particle exits from +y
    inds_plus_y = find(stat_par(:,2) > 1);
    r_par_new(inds_plus_y,2) = mod(r_par_old(inds_plus_y,2), ...
        dom_size(inds_plus_y,2));
elseif ~ isempty(find(stat_par(:,2) < 0,1)) % checking particle exits from -y
    inds_minus_y = find(stat_par(:,2) < 0);
    r_par_new(inds_minus_y,2) = dom_size(inds_minus_y,2) - ...
        mod(abs(r_par_old(inds_minus_y,2)), dom_size(inds_minus_y,2));
end

if ~ isempty(find(stat_par(:,3) > 1,1)) % checking particle exits from +z
    inds_plus_z = find(stat_par(:,3) > 1);
    r_par_new(inds_plus_z,3) = mod(r_par_old(inds_plus_z,3), ...
        dom_size(inds_plus_z,3));
elseif ~ isempty(find(stat_par(:,3) < 0,1)) % checking particle exits from -z
    inds_minus_z = find(stat_par(:,3) < 0);
    r_par_new(inds_minus_z,3) = dom_size(inds_minus_z,3) - ...
        mod(abs(r_par_old(inds_minus_z,3)), dom_size(inds_minus_z,3));
end

end
