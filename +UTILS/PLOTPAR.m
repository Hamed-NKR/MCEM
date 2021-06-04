function [fig_main, varargout] = PLOTPAR(dom_size, par, varargin)
% "PLOTPAR" plots the instantaneous location of primary particles.
% ----------------------------------------------------------------------- %

% Inputs:
    % dom_size: Computational domain dimensions
    % par: The information structure of particles population
    % varargin: An optional varying argument that contains the...
        % ...visibility status of various particle properties:
        % 'equivalent_volumetric_size': 'on'/'off' (or: 'ON/OFF', 'On/Off')
        % 'velocity_vector': 'on'/'off' (or: ~)
        % 'nearest_neighbor': 'on'/'off' (or: ~)
        % 'target_index'(in the case of 'nearest_neighbor' being 'on'):...
            % ...The indices of target particles, an N*1 array of...
            % ...integers (see "COL.NNS" for more info)
        % 'target_coefficient'(again for 'nearest_neighbor' set to...
            % ...'on'): The neighboring limit coefficient (look up...
            % ..."COL.NNS" for details)
        % Note.1: Visibility inputs MUST match the above-mentioned format.
        % Note.2: While 'nearest_neighbor' is 'on'
% ----------------------------------------------------------------------- %

% Outputs
    % fig_main: Figure handle for spatial distribution of particles, and...
        % ...possibly their equivalent size and velocity
    % varargout: Figure handle for the nearest neighbor plots
% ----------------------------------------------------------------------- %        

% Compile primary particles across multiple aggregates.
if isa(par, 'AGG')
    pp = AGG.COMPILEPP(par);
else
    pp = cell2mat(par.pp);
end

% Initializing the main figure handle
hold off
fig_main = gcf;

% Initializing visibility variables
vis_equiv = 0; % Equivalent size visibility status
vis_vel = 0; % Velocity visibility status
vis_nn = 0; % Nearest neighbor visibility status

if nargout > 2
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
    if (mod(nargin,2) ~= 0) || (nargin > 12)
        error('Error: Invalid number of input arguments!')

    elseif numel(unique(varargin_spec)) ~= nargin_spec
        error('Error: Repeating specifications!')

    else
        % Assigning the status variables
        for i = 1 : nargin_spec
            switch varargin_spec{i}
                
                % Checking the status of aggregate equivalents visibility
                case 'equivalent_volumetric_size'
                    vis_equiv = strcmp(varargin{2*i},'on') ||...
                        strcmp(varargin{2*i},'ON') ||...
                        strcmp(varargin{2*i},'On');
                    % Checking for invalid status variables
                    if ~ (vis_equiv || strcmp(varargin{2*i},'off') ||...
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
                    
                % Checking the status of nearest neighbors visibility    
                case 'nearest_neighbor'
                    vis_nn = strcmp(varargin{2*i},'on') ||...
                        strcmp(varargin{2*i},'ON') ||...
                        strcmp(varargin{2*i},'On');
                    if ~ (vis_nn || strcmp(varargin{2*i},'off') ||...
                        strcmp(varargin{2*i},'OFF') ||...
                        strcmp(varargin{2*i},'Off'))
                        error('Error: Invalid input argument number %d \n',...
                            2*i)
                    elseif (~ vis_nn) && (nargout > 1)
                        error('Error: Invalid number of output arguments!')
                    % Initializing the nearest neighbor figure
                    elseif vis_nn
                        ii1 = 2 * find(ismember(varargin_spec,...
                            'target_index'));
                        ind_trg = varargin{ii1};
                        ii2 = 2 * find(ismember(varargin_spec,...
                            'target_coefficient'));
                        coef_trg = varargin{ii2};
                        varargout{1} = figure;  % Initializing the...
                            % ...nearest neighbors figure handle
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



figure(fig_main)

% Setting figure position and background
if ~all(fig_main.Position == [0, 0, 1000, 892.1])
    fig_main.Position = [0, 0, 1000, 892.1];
end
set(fig_main, 'color', 'white');

% Initialing the colors of different plot elements
cm1 = [0.2,0.3,1]; % Color of primaries
cm2 = [0.2,0.8,1]; % Color of spherical aggregate equivalents
cm3 = [1,0.8,0.2]; % Color of velocity vectors

% XY subplot
subplot(2,2,1)
ca = gca;
delete(ca.Children);
% Plotting the primary particles
viscircles([pp(:,3), pp(:,4)], pp(:,2)./2, 'EnhanceVisibility', false,...
    'Color', cm1, 'LineWidth', 0.5); % Plotting primaries
hold on
% Plotting equivalent spheres representing the aggregates
if vis_equiv
    viscircles([par.r(:,1), par.r(:,2)], (par.d)./2, ...
        'EnhanceVisibility', false, 'Color', cm2, 'LineWidth', 2,...
        'LineStyle', '--');
end
% Plotting velocity vectors
if vis_vel
    quiver(par.r(:,1), par.r(:,2), par.v(:,1), par.v(:,2), 'Color', cm3);
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
if vis_equiv
    viscircles([par.r(:,1), par.r(:,3)], (par.d)./2,...
        'EnhanceVisibility', false, 'Color', cm2, 'LineWidth', 2,...
        'LineStyle', '--');
end
if vis_vel
    quiver(par.r(:,1), par.r(:,3), par.v(:,1), par.v(:,3), 'Color', cm3);
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
if vis_equiv
    viscircles([par.r(:,2), par.r(:,3)], (par.d)./2,...
        'EnhanceVisibility', false, 'Color', cm2, 'LineWidth', 2,...
        'LineStyle', '--');
end
if vis_vel
    quiver(par.r(:,2), par.r(:,3), par.v(:,2), par.v(:,3), 'Color', cm3);
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
scatter3(pp(:,3), pp(:,4), pp(:,5), ((pp(:,2))./2).*2e9,...
    cm1,'filled');
hold on
if vis_equiv
    scatter3(par.r(:,1), par.r(:,2), par.r(:,3), ((par.d)./2).*2e9,...
        cm2);
end
if vis_vel
    quiver3(par.r(:,1), par.r(:,2), par.r(:,3), par.v(:,1), par.v(:,2),...
        par.v(:,3), 'Color', cm3);
end
hold off
axis equal
title('xyz view')
xlabel('x (m)')
ylabel('y (m)')
zlabel('z (m)')
xlim([0 dom_size(1)])
ylim([0 dom_size(2)])
zlim([0 dom_size(3)])

% Plotting the nearest neighbors
if vis_nn
    figure(varargout{1})
    varargout{1} = UTILS.PLOTNN(dom_size, par, ind_trg, coef_trg);
end

end

