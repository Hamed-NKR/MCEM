function PLOTPAR(dom_size,par,vis_vel)
% "PLOTPAR" plots the instantaneous location of primary particles.

% dom_size is the computational domain dimensions array.
% pp is the primaries structure.
% vis_vel is a logical variable indicating the visibility of velocity...
% vectors.

hold off;
fig = gcf;
clf(fig);

cm1 = [0.2,0.3,1];  % color of circles
cm2 = [1,0.3,0.2];  % color of velocity vectors

% XY subplot
subplot(2,2,1)
viscircles([par.r(:,1),par.r(:,2)],(par.d)./2, ...
    'EnhanceVisibility',false,'Color',cm1,'LineWidth',0.5);
if vis_vel
    hold on;
    quiver(par.r(:,1),par.r(:,2),par.v(:,1),par.v(:,2),'Color',cm2);
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
viscircles([par.r(:,1),par.r(:,3)],(par.d)./2, ...
    'EnhanceVisibility',false,'Color',cm1,'LineWidth',0.5);
if vis_vel
    hold on;
    quiver(par.r(:,1),par.r(:,3),par.v(:,1),par.v(:,3),'Color',cm2);
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
viscircles([par.r(:,2),par.r(:,3)],(par.d)./2, ...
    'EnhanceVisibility',false,'Color',cm1,'LineWidth',0.5);
if vis_vel
    hold on;
    quiver(par.r(:,2),par.r(:,3),par.v(:,2),par.v(:,3),'Color',cm2);
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
scatter3(par.r(:,1),par.r(:,2),par.r(:,3),((par.d)./2).*2e8,...
    cm1,'filled');
if vis_vel
    hold on
    quiver3(par.r(:,1),par.r(:,2),par.r(:,3),par.v(:,1),par.v(:,2),par.v(:,3),...
        'Color',cm2)
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

