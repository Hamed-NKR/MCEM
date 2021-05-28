function PLOTPAR(dom_size, par, varargin)
% "PLOTPAR" plots the instantaneous location of primary particles.

% Inputs:
    % dom_size: Computational domain dimensions
    % par: The information structure of particles population
    % varargin: An optional varying argument that contains the...
        % ...visibility status of various particle properties:
        % 'equivalent_volumetric_size': 'on'/'off' (or: 'ON/OFF', 'On/Off')
        % 'velocity_vector': 'on'/'off' (or: ~)
        % 'nearest_neighbor': 'on'/'off' (or: ~)
        % Note: Visibility inputs MUST match the above-mentioned format.

hold off;
fig = gcf;
clf(fig);

cm1 = [0.2,0.3,1];  % Color of primaries
cm2 = [0.2,0.8,1];  % Color of equivalent aggregate
cm3 = [1,0.3,0.2];  % Color of velocity vectors

pp = cell2mat(par.pp);

%%% Visibility status
if nargin > 2
    
    vis_ind = 1 : 2 : (nargin -2) ; % Looping index for initialization...
        % ...of visibility status
    nargin_spec = (nargin - 2) / 2;
    varargin_spec = cell(1,nargin_spec);
    for i = 1 : nargin_spec
        varargin_spec{i} = cell2mat(varargin(vis_ind(i)));
    end
    
    if (mod(nargin,2) ~= 0) || (nargin > 8)
        error('Error: Invalid number of input arguments!')

    elseif size(unique(varargin_spec),2) ~= nargin_spec
        error('Error: Repeating elements!')

    else
        % Assigning status variables
        for i = vis_ind
            switch varargin{i}

                case 'equivalent_volumetric_size'
                    vis_equiv = strcmp(varargin{i+1},'on') ||...
                        strcmp(varargin{i+1},'ON') ||...
                        strcmp(varargin{i+1},'On');
                    if ~ (vis_equiv || strcmp(varargin{i+1},'off') ||...
                        strcmp(varargin{i+1},'OFF') ||...
                        strcmp(varargin{i+1},'Off'))
                        error('Error: Invalid input argument number %d \n',...
                            i+1)
                    end

                case 'velocity_vector'
                    vis_vel = strcmp(varargin{i+1},'on') ||...
                        strcmp(varargin{i+1},'ON') ||...
                        strcmp(varargin{i+1},'On');
                    if ~ (vis_vel || strcmp(varargin{i+1},'off') ||...
                        strcmp(varargin{i+1},'OFF') ||...
                        strcmp(varargin{i+1},'Off'))
                        error('Error: Invalid input argument number %d \n',...
                            i+1)
                    end

                case 'nearest_neighbor'
                    vis_nn = strcmp(varargin{i+1},'on') ||...
                        strcmp(varargin{i+1},'ON') ||...
                        strcmp(varargin{i+1},'On');
                    if ~ (vis_nn || strcmp(varargin{i+1},'off') ||...
                        strcmp(varargin{i+1},'OFF') ||...
                        strcmp(varargin{i+1},'Off'))
                        error('Error: Invalid input argument number %d \n',...
                            i+1)
                    end

                otherwise
                    error('Error: Invalid input argument number %d \n', i)

            end

        end

    end
    
end
%%%

% XY subplot
subplot(2,2,1)

viscircles([pp(:,3), pp(:,4)], pp(:,2)./2, 'EnhanceVisibility', false,...
    'Color', cm1, 'LineWidth', 0.5); % Plotting primaries

% Plotting equivalent spheres representing the aggregates
if vis_equiv
    viscircles([par.r(:,1), par.r(:,2)], (par.d)./2, ...
        'EnhanceVisibility', false, 'Color', cm2, 'LineWidth', 2,...
        'LineStyle', '--');
end

% Plotting velocity vectors
if vis_vel
    hold on;
    quiver(par.r(:,1), par.r(:,2), par.v(:,1), par.v(:,2), 'Color', cm3);
    hold off;
end

axis equal
title('xy view')
xlabel('x (m)')
ylabel('y (m)')
xlim([0 dom_size(1)])
ylim([0 dom_size(2)])

% XZ subplot
subplot(2,2,2)

viscircles([pp(:,3), pp(:,5)], pp(:,2)./2, 'EnhanceVisibility', false,...
    'Color', cm1, 'LineWidth' ,0.5);

if vis_equiv
    viscircles([par.r(:,1), par.r(:,3)], (par.d)./2,...
        'EnhanceVisibility', false, 'Color', cm2, 'LineWidth', 2,...
        'LineStyle', '--');
end

if vis_vel
    hold on;
    quiver(par.r(:,1), par.r(:,3), par.v(:,1), par.v(:,3), 'Color', cm3);
    hold off;
end

axis equal
title('xz view')
xlabel('x (m)')
ylabel('z (m)')
xlim([0 dom_size(1)])
ylim([0 dom_size(3)])

% YZ subplot
subplot(2,2,3)

viscircles([pp(:,4), pp(:,5)], pp(:,2)./2, 'EnhanceVisibility', false,...
    'Color', cm1, 'LineWidth' ,0.5);

if vis_equiv
    viscircles([par.r(:,2), par.r(:,3)], (par.d)./2,...
        'EnhanceVisibility', false, 'Color', cm2, 'LineWidth', 2,...
        'LineStyle', '--');
end

if vis_vel
    hold on;
    quiver(par.r(:,2), par.r(:,3), par.v(:,2), par.v(:,3), 'Color', cm3);
    hold off;
end

axis equal
title('yz view')
xlabel('y (m)')
ylabel('z (m)')
xlim([0 dom_size(2)])
ylim([0 dom_size(3)])

% 3D subplot
subplot(2,2,4)

scatter3(pp(:,3), pp(:,4), pp(:,5), ((pp(:,2))./2).*2e9,...
    cm1,'filled');
hold on

if vis_equiv
    scatter3(par.r(:,1), par.r(:,2), par.r(:,3), ((par.d)./2).*2e9,...
        cm2);
end

if vis_vel
    quiver3(par.r(:,1), par.r(:,2), par.r(:,3), par.v(:,1), par.v(:,2),...
        par.v(:,3), 'Color', cm3)
    hold off;
end

axis equal
title('xyz view')
xlabel('x (m)')
ylabel('y (m)')
zlabel('z (m)')
xlim([0 dom_size(1)])
ylim([0 dom_size(2)])
zlim([0 dom_size(3)])

end

