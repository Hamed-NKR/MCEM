function h_pp = PLOTPP(x_pp, y_pp, z_pp, d_pp, fc, ft)
% "PLOTPP" displays the 3d schematic of primary particles in 3d space.
% ----------------------------------------------------------------------- %
% 
% Input:
%     x_pp: x coordinate of primaries 
%     y_pp: y ~
%     z_pp: z ~
%     d_pp: Dimater of primaries
%     fc: Particles color
%     ft: ~ transparency
% ----------------------------------------------------------------------- %
% 
% Output:
%     h_pp: The output figure handle
% ----------------------------------------------------------------------- %

% Setting default colormap as "summer"
if ~exist('fc', 'var'); fc = []; end
if isempty(fc); fc = summer; end


% Setting default marker transparency as 1 (fully opaque)
if ~exist('ft', 'var'); ft = []; end
if isempty(ft); ft = 1; end

n_pp = numel(x_pp); % Number of primary particles

% Plotting aggregates
[X,Y,Z] = sphere(60);

for i = 1 : n_pp
    h_pp = surf(X .* d_pp(i) ./ 2 + x_pp(i),...
        Y .* d_pp(i) ./ 2 + y_pp(i), Z .* d_pp(i) ./ 2 + z_pp(i)); % plot primaries
    
    % set graphics
    h_pp.EdgeColor = 'none';
    h_pp.FaceColor = fc;
    h_pp.FaceAlpha = ft;
    h_pp.FaceLighting = 'gouraud';
    h_pp.AmbientStrength = 0.8;
    h_pp.DiffuseStrength = 0.2;
    h_pp.SpecularStrength = 0.05;
    h_pp.SpecularExponent = 2;
    h_pp.BackFaceLighting = 'lit';
    
    hold on
end

axis equal
grid off
axis off
camlight('right');

if nargout == 0
    clear h_pp;
end

end

