function PLOTPAR(dom_size , par, vis)
% "PLOTPAR" plots the instantaneous location of primary particles.

% Inputs:
    % dom_size: Computational domain dimensions
    % par: The information structure of particles population
    % vis: The visibility status of various particle properties...
        % ...(equivalent aggregate size, velocity vectors, nearest...
        % ...neighbors)

hold off;
fig = gcf;
clf(fig);

cm1 = [0.2,0.3,1];  % Color of primaries
cm2 = [0.2,0.8,1];  % Color of equivalent aggregate
cm3 = [1,0.3,0.2];  % Color of velocity vectors

pp = cell2mat(par.pp);

% XY subplot
subplot(2,2,1)
viscircles([pp(:,3), pp(:,4)], pp(:,2)./2, 'EnhanceVisibility', false,...
    'Color', cm1, 'LineWidth', 0.5); % Plotting primaries
viscircles([par.r(:,1), par.r(:,2)], (par.d)./2, ...
    'EnhanceVisibility', false, 'Color', cm2, 'LineWidth', 2,...
    'LineStyle', '--'); % Plotting equivalent spheres of the aggregates
if vis
    hold on;
    quiver(par.r(:,1), par.r(:,2), par.v(:,1), par.v(:,2), 'Color', cm3);
    % Plotting velocity vectors
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
viscircles([par.r(:,1), par.r(:,3)], (par.d)./2,...
    'EnhanceVisibility', false, 'Color', cm2, 'LineWidth', 2,...
    'LineStyle', '--');
if vis
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
viscircles([par.r(:,2), par.r(:,3)], (par.d)./2,...
    'EnhanceVisibility', false, 'Color', cm2, 'LineWidth', 2,...
    'LineStyle', '--');
if vis
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
scatter3(par.r(:,1), par.r(:,2), par.r(:,3), ((par.d)./2).*2e9,...
    cm2);
if vis
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

