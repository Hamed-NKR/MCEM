
function DYNPLOT(dom_size,pp_d,pp_r)
% This function plots the instantaneous location of a primary particle.

% dom_size is the computational domain dimensions array.
% n_pp is the number of primary particles, and pp_d and pp_r are their...
% size and location.

hold off;
clf;

cm = [0.2,0.3,1];  % color of circles

pp_r = real(pp_r);

% XY subplot
subplot(2,2,1)
viscircles([pp_r(:,1),pp_r(:,2)],pp_d./2, ...
    'EnhanceVisibility',false,'Color',cm,'LineWidth',0.5);
axis equal
title('xy view')
xlabel('x (m)')
ylabel('y (m)')
xlim([0 dom_size(1)])
ylim([0 dom_size(2)])

% XZ subplot
subplot(2,2,2)
viscircles([pp_r(:,1),pp_r(:,3)],pp_d./2, ...
    'EnhanceVisibility',false,'Color',cm,'LineWidth',0.5);
axis equal
title('xz view')
xlabel('x (m)')
ylabel('z (m)')
xlim([0 dom_size(1)])
ylim([0 dom_size(3)])

% YZ subplot
subplot(2,2,3)
viscircles([pp_r(:,2),pp_r(:,3)],pp_d./2, ...
    'EnhanceVisibility',false,'Color',cm,'LineWidth',0.5);
axis equal
title('yz view')
xlabel('y (m)')
ylabel('z (m)')
xlim([0 dom_size(2)])
ylim([0 dom_size(3)])

% 3D subplot
subplot(2,2,4)
scatter3(pp_r(:,1),pp_r(:,2),pp_r(:,3),(pp_d(:)/2)*10^9,...
    cm,'filled');
axis equal
title('xyz view')
xlabel('x (m)')
ylabel('y (m)')
zlabel('z (m)')
xlim([0 dom_size(1)])
ylim([0 dom_size(2)])
zlim([0 dom_size(3)])

end

