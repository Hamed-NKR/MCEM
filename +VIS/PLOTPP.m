
function PLOTPP(dom_size,pp,vis_vel)
% This function plots the instantaneous location of primary particles.

% dom_size is the computational domain dimensions array.
% pp is the primaries structure.
% vis_vel is a logical variable indicating the visibility of velocity...
% vectors.

hold off;
clf;

cm1 = [0.2,0.3,1];  % color of circles
cm2 = [1,0.3,0.2];  % color of velocity vectors

pp_r = real(pp.r);

% XY subplot
subplot(2,2,1)
viscircles([pp_r(:,1),pp_r(:,2)],(pp.d)./2, ...
    'EnhanceVisibility',false,'Color',cm1,'LineWidth',0.5);
if vis_vel
    hold on;
    quiver(pp.r(:,1),pp.r(:,2),pp.v(:,1),pp.v(:,2));
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
viscircles([pp_r(:,1),pp_r(:,3)],(pp.d)./2, ...
    'EnhanceVisibility',false,'Color',cm1,'LineWidth',0.5);
if vis_vel
    hold on;
    quiver(pp.r(:,1),pp.r(:,3),pp.v(:,1),pp.v(:,3));
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
viscircles([pp_r(:,2),pp_r(:,3)],(pp.d)./2, ...
    'EnhanceVisibility',false,'Color',cm1,'LineWidth',0.5);
if vis_vel
    hold on;
    quiver(pp.r(:,2),pp.r(:,3),pp.v(:,2),pp.v(:,3));
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
scatter3(pp_r(:,1),pp_r(:,2),pp_r(:,3),((pp.d)./2)*10^9,...
    cm1,'filled');
if vis_vel
    hold on
    quiver3(pp.r(:,1),pp.r(:,2),pp.r(:,3),pp.v(:,1),pp.v(:,2),pp.v(:,3))
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

