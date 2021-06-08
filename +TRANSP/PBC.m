function par = PBC(dom_size, par)
% "PBC" applies periodic boundary condition to the particle movements.
% ----------------------------------------------------------------------- %

% Inputs/Outputs:
    % dom_size: Computational domain size
    % par: Particle information structure
% ----------------------------------------------------------------------- %

if isa(par, 'AGG')
    n_par = size(par, 1);
else
    n_par = size(par.n, 1); % Total number of particles
end

r = cat(1, par.r);

dom_size = repmat(dom_size', n_par, 1);
stat_par = r ./ dom_size; % This status array shows whether the...
    % ...particles are inside the domain.

% Finding the new positions for the case of domain escape
par_r_new = r; % Temporary array for new positions
if ~ isempty(find(stat_par(:,1) > 1,1)) % checking particle exits from +x
    inds_plus_x = find(stat_par(:,1) > 1); % finding the out of domain...
        % ...particles
    par_r_new(inds_plus_x,1) = mod(r(inds_plus_x,1),...
        dom_size(inds_plus_x,1)); % Replacing the outgoing particles...
            % ...with new ones coming from the opposite side.
elseif ~ isempty(find(stat_par(:,1) < 0,1)) % checking particle exits...
        % ...from -x
    inds_minus_x = find(stat_par(:,1) < 0);
    par_r_new(inds_minus_x,1) = dom_size(inds_minus_x,1) -...
        mod(abs(r(inds_minus_x,1)), dom_size(inds_minus_x,1));
end

if ~ isempty(find(stat_par(:,2) > 1,1)) % checking particle exits from +y
    inds_plus_y = find(stat_par(:,2) > 1);
    par_r_new(inds_plus_y,2) = mod(r(inds_plus_y,2),...
        dom_size(inds_plus_y,2));
elseif ~ isempty(find(stat_par(:,2) < 0,1)) % checking particle exits...
        % ...from -y
    inds_minus_y = find(stat_par(:,2) < 0);
    par_r_new(inds_minus_y,2) = dom_size(inds_minus_y,2) -...
        mod(abs(r(inds_minus_y,2)), dom_size(inds_minus_y,2));
end

if ~ isempty(find(stat_par(:,3) > 1,1)) % checking particle exits from +z
    inds_plus_z = find(stat_par(:,3) > 1);
    par_r_new(inds_plus_z,3) = mod(r(inds_plus_z,3),...
        dom_size(inds_plus_z,3));
elseif ~ isempty(find(stat_par(:,3) < 0,1)) % checking particle exits...
        % ...from -z
    inds_minus_z = find(stat_par(:,3) < 0);
    par_r_new(inds_minus_z,3) = dom_size(inds_minus_z,3) -...
        mod(abs(r(inds_minus_z,3)), dom_size(inds_minus_z,3));
end

dr = par_r_new - r;

if isa(par, 'AGG')
    par = par.TRANSLATE(dr);

else    
    % Updating the primary particle locations
    pp = cell2mat(par.pp);
    pp(:,3:5) = pp(:,3:5) + repelem(dr, par.n, 1);
    par.pp = mat2cell(pp, par.n);

    % Updating the particle center of mass
    par.r = par_r_new;
end

end
