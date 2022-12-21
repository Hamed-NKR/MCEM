function h_cc = COLORCODE(pars, dom_size, cm, q)
% "COLORCODE" visualizes hybrid aggregates with nonunform primaries...
%   ...using different colors for different regions of the near-size
%   ...primaries.
% ----------------------------------------------------------------------- %
% 
% Inputs:
%     pars: Particle information structure/class
%     dom_size: Computational domain dimensions
%     cm: Plot colormap
%     q: Quality of render.
% ----------------------------------------------------------------------- %
% 
% Output:
%     h_cc: The output figure handle
% ----------------------------------------------------------------------- %

warning('off')

% Setting default domain size as empty
if ~exist('dom_size', 'var'); dom_size = []; end

% Setting default colormap as "autumn"
if ~exist('cm', 'var'); cm = []; end
if isempty(cm); cm = autumn; end

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
h_cc = gcf;
figure(h_cc);
if ~all(h_cc.Position == [0, 0, 1000, 892.1])
    h_cc.Position = [0, 0, 1000, 892.1]; % Setting position
end
set(h_cc, 'color', 'white');

if isa(pars, 'AGG')
    pp = AGG.COMPILEPP(pars);
else
    pp = cell2mat(pars.pp);
end
pp_d = unique(pp(:,2));
ii = unique(round(10 + (length(cm) - 10) .*...
    (0 : 1 / (length(pp_d) - 1) : 1)'));
cl = cm(ii,:);

% Plotting aggregates
[X,Y,Z] = sphere(q);

for i = 1 : n_par    
    for j = 1 : n_pp(i)
        
        if isa(pars, 'AGG')
            h_cc = surf(X .* pars(i).pp.d(j) ./ 2 + pars(i).pp.r(j,1),...
                Y .* pars(i).pp.d(j) ./ 2 + pars(i).pp.r(j,2),...
                Z .* pars(i).pp.d(j) ./ 2 + pars(i).pp.r(j,3));
            i_cl = round((find(pp_d == pars(i).pp.d(j),1) /...
                 numel(pp_d)) * length(cl));
            h_cc.FaceColor = cl(i_cl,:);
        else
            h_cc = surf(X .* pars.pp{i}(j,2) ./ 2 + pars.pp{i}(j,3),...
                Y .* pars.pp{i}(j,2) ./ 2 + pars.pp{i}(j,4),...
                Z .* pars.pp{i}(j,2) ./ 2 + pars.pp{i}(j,5));
            i_cl = find(pp_d == pars.pp{i}(j,2));
            i_cl = round(i_cl(randperm(numel(i_cl),1)) / numel(pp_d)) *...
                length(cl);
            h_cc.FaceColor = cl(i_cl,:);
        end
        
        lightangle(-45,30)
        h_cc.EdgeColor = 'none';
        h_cc.FaceLighting = 'gouraud';
        h_cc.AmbientStrength = 0.7;
        h_cc.DiffuseStrength = 0.1;
        h_cc.SpecularStrength = 0.01;
        h_cc.SpecularExponent = 2;
        h_cc.BackFaceLighting = 'lit';
        
        
        hold on
        UTILS.TEXTBAR([sum(n_pp(1 : i-1)) + j, n_tot]);
    end
end

disp(' ');
        
% Formatting plot
disp('Formatting plot ...');

view(3);                                  % view to start from
c = camlight('headlight');                % add light
set(c,'style','infinite');                % set style of light

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

