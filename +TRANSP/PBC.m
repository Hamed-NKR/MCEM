function pars = PBC(dom_size, pars)
% "PBC" applies periodic boundary condition to the particle movements.
% ----------------------------------------------------------------------- %
% 
% Inputs/Outputs:
%     dom_size: Computational domain size
%     pars: The particle information structure/class
% ----------------------------------------------------------------------- %

% Total number of (independent) particles
if isa(pars, 'AGG')
    n_par = length(pars);
else
    n_par = size(pars.n, 1);
end

% Compiling/copying properties locally
r = cat(1, pars.r);

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

dr = par_r_new - r; % Translation vectors

% Updating the locations for the structure/class of aggregates
if isa(pars, 'AGG')
    pars = pars.TRANSLATE(dr);

else
    [pars.pp, pars.r] = PAR.TRANSLATE(pars.pp, pars.r, pars.n, dr);
end

end
