function [r_cm, varargout] = EQUIV(pp, n_pp)
% "EQUIV" computes the center of mass and equivalent volumetric size of...
%     ...an ensemble of spherical particles.
% ----------------------------------------------------------------------- %
% 
% Inputs:
%     pp: Primary particle information cell array
%     n_pp: Number distribution of primaries
% ----------------------------------------------------------------------- %
% 
% Outputs:
%     r_cm: Aggregate center of mass position
%     varargout{1}: Aggregate equivalent volumetric diameter
% ----------------------------------------------------------------------- %

if nargout > 2
    error('Error: Invalid number of output arguments!') % Checking for...
        % ...redundant output arguments
end

pp = cell2mat(pp);
vol_pp = pi .* (pp(:,2).^3) ./ 6; % Primary particle volumes
r_temp = (vol_pp .* pp(:,3:5));
vol_pp = mat2cell(vol_pp, n_pp);
r_temp = mat2cell(r_temp, n_pp);
n_agg = size(n_pp,1); % Total number of aggregates
r_cm = zeros(n_agg,3);
if nargout > 1
    varargout{1} = zeros(n_agg,1);
end

for i = 1 : n_agg
    r_cm(i,:) = sum(r_temp{i}) ./ sum(vol_pp{i}); % Calculating the...
        % ...center of mass
    if nargout > 1
        varargout{1}(i) = nthroot(6 .* sum(vol_pp{i}) ./ pi, 3);
            % Calculating the equivalent diameter
    end
end

end
