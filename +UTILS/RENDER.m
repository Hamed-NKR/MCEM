function h = RENDER(dom_size, pps, n_pps, cm)
% "RENDER" depicts the posture of an aggregate population in space.
% 
% NOTE: Originally developed by Tim Sipkens, 2021.
% ----------------------------------------------------------------------- %
% 
% Input:
%     dom_size: Computational domain dimensions
%     pps: A cell array of primary particles info
%     n_pps: Number distribution of primaries
%     cm: Plot colormap
% ----------------------------------------------------------------------- %
% 
% Output:
%     fig: The output figure handle
% ----------------------------------------------------------------------- %

% Setting default colormap as "summer"
if ~exist('cm', 'var'); cm = []; end
if isempty(cm); cm = summer; end

n_par = numel(pps); % Number of independent particles
n_tot = sum(n_pps); % Total number of primaries

% Displaying overal progress
disp('Rendering:');
UTILS.TEXTBAR([0, n_tot]);

% Initializing the figure handle
hold off
h = gcf;
figure(h);
if ~all(h.Position == [0, 0, 1000, 892.1])
    h.Position = [0, 0, 1000, 892.1]; % Setting position
end
set(h, 'color', 'white');
colormap(cm); % Setting colormap

% Plotting aggregates
[X,Y,Z] = sphere(60);
for i = 1 : n_par
        
    for j = 1 : n_pps(i)
    
        h = surf(X .* pps{i}(j,2) ./ 2 + pps{i}(j,3), ...
            Y .* pps{i}(j,2) ./ 2 + pps{i}(j,4), ...
            Z .* pps{i}(j,2) ./ 2 + pps{i}(j,5));
        
        hold on
        UTILS.TEXTBAR([sum(n_pps(1 : i-1)) + j, n_tot]);
    end
end

lightangle(-45,30)
h.FaceLighting = 'gouraud';
h.AmbientStrength = 0.8;
h.DiffuseStrength = 0.2;
h.SpecularStrength = 0.05;
h.SpecularExponent = 2;
h.BackFaceLighting = 'lit';

disp(' ');
        
% Formatting plot
disp('Formatting plot ...');
camlight('right');
shading interp;
view([-37, 20]);
hold off
axis equal;
xlabel('x (m)')
ylabel('y (m)')
zlabel('z (m)')
xlim([0 dom_size(1)])
ylim([0 dom_size(2)])
zlim([0 dom_size(3)])

disp('DONE.');
disp(' ');

if nargout == 0
    clear h;
end

end

