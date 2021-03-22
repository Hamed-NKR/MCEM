function DYNPLOT(dom_size,n_pp,pp_d,pp_r)
% This function plots the instantaneous location of a primary particle.

% dom_size is the computational domain dimensions array.
% n_pp is the number of primary particles, and pp_d and pp_r are their...
% size and location.

subplot(2,2,1)
for j = 1:n_pp
    VIS.CIRCLE(pp_r(j,1),pp_r(j,2),pp_d(j)/2);
    hold on
end
axis equal
title('xy view')
xlabel('x (m)')
ylabel('y (m)')
xlim([0 dom_size(1)])
ylim([0 dom_size(2)])
hold off

subplot(2,2,2)
for j = 1:n_pp
    VIS.CIRCLE(pp_r(j,1),pp_r(j,3),pp_d(j)/2);
    hold on
end
axis equal
title('xz view')
xlabel('x (m)')
ylabel('z (m)')
xlim([0 dom_size(1)])
ylim([0 dom_size(3)])
hold off

subplot(2,2,3)
for j = 1:n_pp
    VIS.CIRCLE(pp_r(j,2),pp_r(j,3),pp_d(j)/2);
    hold on
end
axis equal
title('yz view')
xlabel('y (m)')
ylabel('z (m)')
xlim([0 dom_size(2)])
ylim([0 dom_size(3)])
hold off

subplot(2,2,4)
% for j = 1:n_pp
%     VIS.SPHERE(pp_r(j,1),pp_r(j,2),pp_r(j,3),pp_d(j)/2);
%     hold on
% end
scatter3(pp_r(:,1),pp_r(:,2),pp_r(:,3),(pp_d(:)/2)*10^9);
hold on
axis equal
title('xyz view')
xlabel('x (m)')
ylabel('y (m)')
zlabel('z (m)')
xlim([0 dom_size(1)])
ylim([0 dom_size(2)])
zlim([0 dom_size(3)])
hold off

end

