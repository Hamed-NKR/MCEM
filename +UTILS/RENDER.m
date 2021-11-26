function h_rend = RENDER(pars, dom_size, cm, q)
% "RENDER" depicts the posture of an aggregate population in space.
% 
% Original author: Timothy Sipkens, 05-2021
% Revised by: Hamed Nikookar, 06-2021
% ----------------------------------------------------------------------- %
% 
% Input:
%     pars: Particle information structure/class
%     dom_size: Computational domain dimensions
%     cm: Plot colormap
%     q: Quality of render.
% ----------------------------------------------------------------------- %
% 
% Output:
%     h_rend: The output figure handle
% ----------------------------------------------------------------------- %

warning('off')

% Setting default domain size as empty
if ~exist('dom_size', 'var'); dom_size = []; end

% Setting default colormap as "summer"
if ~exist('cm', 'var'); cm = []; end
if isempty(cm); cm = summer; end

% Setting default quality to 60 points in the sphere.
if ~exist('q', 'var'); q = []; end
if isempty(q); q = 60; end

% Compiling primary particles across multiple aggregates
if isa(pars, 'AGG')
    n_par = length(pars); % Number of independent particles
else
    n_par = length(pars.n);
end

n_pp = cat(1, pars.n); % Compiling/copying number distribution of primaries
n_tot = sum(n_pp); % Total number of primaries

% Displaying overal progress
disp(' ')
disp('Rendering:')
UTILS.TEXTBAR([0, n_tot]);

% Initializing the figure handle
hold off
h_rend = gcf;
figure(h_rend);
if ~all(h_rend.Position == [0, 0, 1000, 892.1])
    h_rend.Position = [0, 0, 1000, 892.1]; % Setting position
end
set(h_rend, 'color', 'white');
colormap(cm); % Setting colormap

% Plotting aggregates
[X,Y,Z] = sphere(q);

for i = 1 : n_par    
    for j = 1 : n_pp(i)
        
        if isa(pars, 'AGG')
            h_rend = surf(X .* pars(i).pp.d(j) ./ 2 + pars(i).pp.r(j,1),...
                Y .* pars(i).pp.d(j) ./ 2 + pars(i).pp.r(j,2),...
                Z .* pars(i).pp.d(j) ./ 2 + pars(i).pp.r(j,3));
        
        else
            h_rend = surf(X .* pars.pp{i}(j,2) ./ 2 + pars.pp{i}(j,3),...
                Y .* pars.pp{i}(j,2) ./ 2 + pars.pp{i}(j,4),...
                Z .* pars.pp{i}(j,2) ./ 2 + pars.pp{i}(j,5));
        end
        
        lightangle(-45,30)
        h_rend.FaceLighting = 'gouraud';
        h_rend.AmbientStrength = 0.8;
        h_rend.DiffuseStrength = 0.2;
        h_rend.SpecularStrength = 0.05;
        h_rend.SpecularExponent = 2;
        h_rend.BackFaceLighting = 'lit';
        
        hold on
        UTILS.TEXTBAR([sum(n_pp(1 : i-1)) + j, n_tot]);
    end
end

disp(' ');
        
% Formatting plot
disp('Formatting plot ...');
camlight('right');
shading interp;
view([-37, 20]);
hold off
axis equal;
if isempty(dom_size)
    axis('off')
else
    xlabel('x (m)')
    ylabel('y (m)')
    zlabel('z (m)')
    xlim([0 dom_size(1)])
    ylim([0 dom_size(2)])
    zlim([0 dom_size(3)])
    grid off
end

disp('DONE.');
disp(' ');

if nargout == 0
    clear h_rend; % Deleting figure handle if not requested as an output
end

end

