function STCPLOT(dom_size,pp)
% This function plots the instantaneous location and velocity of...
% the primaries.

% dom_size is the computational domain dimensions array.
% pp is the primary particles' structure.

n_pp = size(pp.ind,1);

subplot(2,2,1)
for j = 1:n_pp
    VIS.CIRCLE(pp.r(j,1),pp.r(j,2),pp.d(j)/2);
    hold on
end
quiver(pp.r(:,1),pp.r(:,2),pp.v(:,1),pp.v(:,2));
axis equal
title('xy view')
xlabel('x (m)')
ylabel('y (m)')
xlim([0 dom_size(1)])
ylim([0 dom_size(2)])

subplot(2,2,2)
for j = 1:n_pp
    VIS.CIRCLE(pp.r(j,1),pp.r(j,3),pp.d(j)/2);
    hold on
end
quiver(pp.r(:,1),pp.r(:,3),pp.v(:,1),pp.v(:,3));
axis equal
title('xz view')
xlabel('x (m)')
ylabel('z (m)')
xlim([0 dom_size(1)])
ylim([0 dom_size(3)])

subplot(2,2,3)
for j = 1:n_pp
    VIS.CIRCLE(pp.r(j,2),pp.r(j,3),pp.d(j)/2);
    hold on
end
quiver(pp.r(:,2),pp.r(:,3),pp.v(:,2),pp.v(:,3));
axis equal
title('yz view')
xlabel('y (m)')
ylabel('z (m)')
xlim([0 dom_size(2)])
ylim([0 dom_size(3)])

subplot(2,2,4)
% for j = 1:n_pp
%     VIS.SPHERE(pp.r(j,1),pp.r(j,2),pp.r(j,3),pp.d(j)/2);
%     hold on
% end
scatter3(pp.r(:,1),pp.r(:,2),pp.r(:,3),(pp.d(:)/2)*10^9);
hold on
quiver3(pp.r(:,1),pp.r(:,2),pp.r(:,3),pp.v(:,1),pp.v(:,2),pp.v(:,3))
axis equal
title('xyz view')
xlabel('x (m)')
ylabel('y (m)')
zlabel('z (m)')
xlim([0 dom_size(1)])
ylim([0 dom_size(2)])
zlim([0 dom_size(3)])

end

