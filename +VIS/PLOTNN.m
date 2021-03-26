function PLOTNN(dom_size,par,ind_trg)
% This function plots the nearest neighbors of a desired particles.

% Inputs are computational domain size, particle population structure,...
% ...and target particle index.

hold off;
fig = gcf;
clf(fig);

cm1 = [0.2,1,0.3];  % Color of target particle
cm2 = [1,0.3,0.2];  % color of nearest neighbors
cm3 = [0.3,0.2,1];  % Color of far particles

par_r = real(par.r);
n_par = size(par_r,1); % Getting number of the particles 

ind_nn = COL.NNS(par_r,par.d,ind_trg); % Getting the nearest neighbors of...
% ...the target particle
ind_rest = (1:n_par)';
ind_rest([ind_trg;ind_nn]) = [];

% XY subplot
subplot(2,2,1)
viscircles([par_r(ind_trg,1),par_r(ind_trg,2)],(par.d(ind_trg))./2, ...
    'EnhanceVisibility',false,'Color',cm1,'LineWidth',0.5);
hold on
viscircles([par_r(ind_nn,1),par_r(ind_nn,2)],(par.d(ind_nn))./2, ...
    'EnhanceVisibility',false,'Color',cm2,'LineWidth',0.5);
viscircles([par_r(ind_rest,1),par_r(ind_rest,2)],(par.d(ind_rest))./2, ...
    'EnhanceVisibility',false,'Color',cm3,'LineWidth',0.5);
hold off
axis equal
title('xy view')
xlabel('x (m)')
ylabel('y (m)')
xlim([0 dom_size(1)])
ylim([0 dom_size(2)])

% XZ subplot
subplot(2,2,2)
viscircles([par_r(ind_trg,1),par_r(ind_trg,3)],(par.d(ind_trg))./2, ...
    'EnhanceVisibility',false,'Color',cm1,'LineWidth',0.5);
hold on
viscircles([par_r(ind_nn,1),par_r(ind_nn,3)],(par.d(ind_nn))./2, ...
    'EnhanceVisibility',false,'Color',cm2,'LineWidth',0.5);
viscircles([par_r(ind_rest,1),par_r(ind_rest,3)],(par.d(ind_rest))./2, ...
    'EnhanceVisibility',false,'Color',cm3,'LineWidth',0.5);
hold off
axis equal
title('xz view')
xlabel('x (m)')
ylabel('z (m)')
xlim([0 dom_size(1)])
ylim([0 dom_size(3)])

% YZ subplot
subplot(2,2,3)
viscircles([par_r(ind_trg,2),par_r(ind_trg,3)],(par.d(ind_trg))./2, ...
    'EnhanceVisibility',false,'Color',cm1,'LineWidth',0.5);
hold on
viscircles([par_r(ind_nn,2),par_r(ind_nn,3)],(par.d(ind_nn))./2, ...
    'EnhanceVisibility',false,'Color',cm2,'LineWidth',0.5);
viscircles([par_r(ind_rest,2),par_r(ind_rest,3)],(par.d(ind_rest))./2, ...
    'EnhanceVisibility',false,'Color',cm3,'LineWidth',0.5);
hold off
axis equal
title('yz view')
xlabel('y (m)')
ylabel('z (m)')
xlim([0 dom_size(2)])
ylim([0 dom_size(3)])

% 3D subplot
subplot(2,2,4)
scatter3(par_r(ind_trg,1),par_r(ind_trg,2),par_r(ind_trg,3),...
    ((par.d(ind_trg))./2).*2e9,cm1,'filled');
hold on
scatter3(par_r(ind_nn,1),par_r(ind_nn,2),par_r(ind_nn,3),...
    ((par.d(ind_nn))./2).*2e9,cm2,'filled');
scatter3(par_r(ind_rest,1),par_r(ind_rest,2),par_r(ind_rest,3),...
    ((par.d(ind_rest))./2).*2e9,cm3,'filled');
hold off
axis equal
title('xyz view')
xlabel('x (m)')
ylabel('y (m)')
zlabel('z (m)')
xlim([0 dom_size(1)])
ylim([0 dom_size(2)])
zlim([0 dom_size(3)])

end

