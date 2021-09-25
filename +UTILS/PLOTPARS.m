function h_pars = PLOTPARS(pars, dom_size, varargin)
% "PLOTPAR" plots the instantaneous location of primary particles.
% ----------------------------------------------------------------------- %
%
% Inputs:
%     dom_size: Computational domain dimensions
%     pars: The particle information structure/class
%     varargin: An optional varying argument that contains the...
%         ...visibility status of various particle properties:
%         (a) 'equivalent_size': 'on'/'off' (or: 'ON/OFF', 'On/Off')
%         (b) 'velocity_vector': 'on'/'off' (or: ~)
%
% Note 1: Visibility inputs MUST match the above-mentioned format.
% Note 2: Equivalent size displayed here is the maximum extent of the
%   ...aggregates.
% ----------------------------------------------------------------------- %
%
% Outputs:
%     h_main: Figure handle for spatial distribution of particles, and...
%         ...possibly their equivalent size and velocity
% ----------------------------------------------------------------------- %        

% Initializing the main figure handle
hold off
h_pars = gcf;

% Initializing visibility variables
vis_equi = 0; % Equivalent size visibility status
vis_vel = 0; % Velocity ~
vis_rend = 0; % Render ~

if nargout > 1
    error('Error: Invalid number of output arguments!') % Checking for...
        % ...redundant output arguments
end

%%% initialization of visibility status parameters
if nargin > 2
    
    vis_ind = 1 : 2 : (nargin -2) ; % Looping index for  visibility...
        % ...checking 
    nargin_spec = (nargin - 2) / 2; % Number of visibility...
        % ...specification variables
    varargin_spec = cell(1,nargin_spec); % Visibility specification...
        % ...variables
    for i = 1 : nargin_spec
        varargin_spec{i} = cell2mat(varargin(vis_ind(i)));
    end
    
    % Checking for redundant, insufficient, or repeating input arguments
    if (mod(nargin,2) ~= 0) || (nargin > 8)
        error('Error: Invalid number of input arguments!')

    elseif numel(unique(varargin_spec)) ~= nargin_spec
        error('Error: Repeating specifications!')

    else
        % Assigning the status variables
        for i = 1 : nargin_spec
            switch varargin_spec{i}
                
                % Checking the status of aggregate equivalents visibility
                case 'equivalent_size'
                    vis_equi = strcmp(varargin{2*i},'on') ||...
                        strcmp(varargin{2*i},'ON') ||...
                        strcmp(varargin{2*i},'On');
                    % Checking for invalid status variables
                    if ~ (vis_equi || strcmp(varargin{2*i},'off') ||...
                        strcmp(varargin{2*i},'OFF') ||...
                        strcmp(varargin{2*i},'Off'))
                        error('Error: Invalid input argument number %d \n',...
                            2*i)
                    end
                    
                % Checking the status of particle velocities visibility    
                case 'velocity_vector'
                    vis_vel = strcmp(varargin{2*i},'on') ||...
                        strcmp(varargin{2*i},'ON') ||...
                        strcmp(varargin{2*i},'On');
                    if ~ (vis_vel || strcmp(varargin{2*i},'off') ||...
                        strcmp(varargin{2*i},'OFF') ||...
                        strcmp(varargin{2*i},'Off'))
                        error('Error: Invalid input argument number %d \n',...
                            2*i)
                    end
                    
                % Checking the rendering status
                case 'render'
                    vis_rend = strcmp(varargin{2*i},'on') ||...
                        strcmp(varargin{2*i},'ON') ||...
                        strcmp(varargin{2*i},'On');
                    % Checking for invalid status variables
                    if ~ (vis_rend || strcmp(varargin{2*i},'off') ||...
                        strcmp(varargin{2*i},'OFF') ||...
                        strcmp(varargin{2*i},'Off'))
                        error('Error: Invalid input argument number %d \n',...
                            2*i)
                    end
                    
                otherwise
                    if ~ (strcmp(varargin_spec{i},'target_index') ||...
                            strcmp(varargin_spec{i},'target_coefficient'))
                        error('Error: Invalid input argument type; Arg. no.: %d \n'...
                            , 2*i + 1 )
                    end

            end

        end

    end
    
end
%%%

figure(h_pars);

% Setting figure position and background
if ~all(h_pars.Position == [0, 0, 1000, 892.1])
    h_pars.Position = [0, 0, 1000, 892.1];
end
set(h_pars, 'color', 'white');

% Initialing the colors of different plot elements
cm1 = [0.15, 0.6, 0.4]; % Color of primaries
cm2 = [0.5, 0.75, 0.4]; % Color of spherical aggregate equivalents
cm3 = [0.9, 0.95, 0.4]; % Color of velocity vectors

% Compiling primary particles across multiple aggregates
if isa(pars, 'AGG')
    pp = AGG.COMPILEPP(pars);
else
    pp = cell2mat(pars.pp);
end

% Concatinating the aggregates global info
if vis_equi
    d = cat(1, pars.dmax);
    r = cat(1, pars.r);
end

if vis_vel
    if ~exist('r', 'var'); r = cat(1, pars.r); end
    v = cat(1, pars.v);
end

% XY subplot
subplot(2,2,1)
ca = gca;
delete(ca.Children);
% Plotting the primary particles
viscircles([pp(:,3), pp(:,4)], pp(:,2)./2, 'EnhanceVisibility', false,...
    'Color', cm1, 'LineWidth', 0.5); % Plotting primaries
hold on
% Plotting equivalent spheres representing the aggregates
if vis_equi
    viscircles([r(:,1), r(:,2)], d ./ 2, 'EnhanceVisibility', false,...
        'Color', cm2, 'LineWidth', 2, 'LineStyle', '--');
end
% Plotting velocity vectors
if vis_vel
    quiver(r(:,1), r(:,2), v(:,1), v(:,2), 'Color', cm3);
end
% Setting different plot graphics
hold off;
axis equal
title('xy view')
xlabel('x (m)')
ylabel('y (m)')
xlim([0 dom_size(1)])
ylim([0 dom_size(2)])

% XZ subplot
subplot(2,2,2)
ca = gca;
delete(ca.Children);
viscircles([pp(:,3), pp(:,5)], pp(:,2)./2, 'EnhanceVisibility', false,...
    'Color', cm1, 'LineWidth' ,0.5);
hold on
if vis_equi
    viscircles([r(:,1), r(:,3)], d ./ 2, 'EnhanceVisibility', false,...
        'Color', cm2, 'LineWidth', 2, 'LineStyle', '--');
end
if vis_vel
    quiver(r(:,1), r(:,3), v(:,1), v(:,3), 'Color', cm3);
end
hold off
axis equal
title('xz view')
xlabel('x (m)')
ylabel('z (m)')
xlim([0 dom_size(1)])
ylim([0 dom_size(3)])

% YZ subplot
subplot(2,2,3)
ca = gca;
delete(ca.Children);
hold on
viscircles([pp(:,4), pp(:,5)], pp(:,2)./2, 'EnhanceVisibility', false,...
    'Color', cm1, 'LineWidth' ,0.5);
if vis_equi
    viscircles([r(:,2), r(:,3)], d ./ 2, 'EnhanceVisibility', false,...
        'Color', cm2, 'LineWidth', 2, 'LineStyle', '--');
end
if vis_vel
    quiver(r(:,2), r(:,3), v(:,2), v(:,3), 'Color', cm3);
end
hold off
axis equal
title('yz view')
xlabel('y (m)')
ylabel('z (m)')
xlim([0 dom_size(2)])
ylim([0 dom_size(3)])

% 3D subplot
subplot(2,2,4)
ca = gca;
delete(ca.Children);
if vis_rend
    UTILS.RENDER(pars, dom_size);
else
    scatter3(pp(:,3), pp(:,4), pp(:,5), ((pp(:,2)) ./ 2) .* 2e9, cm1,...
        'filled');
end
hold on
if vis_equi
    if vis_rend
        UTILS.PLOTPP(r(:,1), r(:,2), r(:,3), d, cm2, 0.3);
    else
        scatter3(r(:,1), r(:,2), r(:,3), (d ./ 2) .* 2e9, cm2);
    end
end
if vis_vel
    quiver3(r(:,1), r(:,2), r(:,3), v(:,1), v(:,2), v(:,3), 'Color', cm3);
end
hold off
axis equal
grid off
title('xyz view')
xlabel('x (m)')
ylabel('y (m)')
zlabel('z (m)')
xlim([0 dom_size(1)])
ylim([0 dom_size(2)])
zlim([0 dom_size(3)])

if nargout == 0
    clear h_pars;  % Deleting figure handle if not requested as an output
end

end

