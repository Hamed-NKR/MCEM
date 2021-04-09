function PLOTNN(dom_size,par,ind_trg,coef_trg)
% This function plots the nearest neighbors of a desired particles.

% Inputs are computational domain size, particle population structure,...
% ...and target particle index.

hold off;
fig = gcf;
clf(fig);

cm1 = [0.2,1,0.3];  % Color of target particle
cm2 = [1,0.3,0.2];  % color of nearest neighbors
cm3 = [0.3,0.2,1];  % Color of far particles

n_par = size(par.r,1); % Getting number of the particles 

ind_nn = COL.NNS(par.r,par.d,ind_trg,coef_trg); % Getting the nearest...
% ...neighbors of the target particle
ind_rest = (1:n_par)';
ind_rest([ind_trg;ind_nn]) = [];

tiledlayout(2,2);

% XY tile
nexttile
p1_xy = viscircles([par.r(ind_trg,1),par.r(ind_trg,2)],(par.d(ind_trg))./2, ...
    'EnhanceVisibility',false,'Color',cm1,'LineWidth',0.5);
hold on
p2_xy = viscircles([par.r(ind_trg,1),par.r(ind_trg,2)],...
    (coef_trg.*(par.d(ind_trg))./2)+(max(par.d)/2),...
    'EnhanceVisibility',false,'Color','k','LineStyle','--','LineWidth',0.5);
p3_xy = viscircles([par.r(ind_nn,1),par.r(ind_nn,2)],(par.d(ind_nn))./2, ...
    'EnhanceVisibility',false,'Color',cm2,'LineWidth',0.5);
p4_xy = viscircles([par.r(ind_rest,1),par.r(ind_rest,2)],(par.d(ind_rest))./2, ...
    'EnhanceVisibility',false,'Color',cm3,'LineWidth',0.5);
hold off
axis equal
title('xy view')
xlabel('x (m)')
ylabel('y (m)')
lgd1 = legend([p1_xy p2_xy p3_xy p4_xy],{'Target','Limit','Neighbor','Rest'});
lgd1.Layout.Tile = 'east';
xlim([0 dom_size(1)])
ylim([0 dom_size(2)])

% XZ tile
nexttile
viscircles([par.r(ind_trg,1),par.r(ind_trg,3)],(par.d(ind_trg))./2, ...
    'EnhanceVisibility',false,'Color',cm1,'LineWidth',0.5);
hold on
viscircles([par.r(ind_trg,1),par.r(ind_trg,3)],...
    (coef_trg.*(par.d(ind_trg))./2)+(max(par.d)/2),...
    'EnhanceVisibility',false,'Color','k','LineStyle','--','LineWidth',0.5);
viscircles([par.r(ind_nn,1),par.r(ind_nn,3)],(par.d(ind_nn))./2, ...
    'EnhanceVisibility',false,'Color',cm2,'LineWidth',0.5);
viscircles([par.r(ind_rest,1),par.r(ind_rest,3)],(par.d(ind_rest))./2, ...
    'EnhanceVisibility',false,'Color',cm3,'LineWidth',0.5);
hold off
axis equal
title('xz view')
xlabel('x (m)')
ylabel('z (m)')
xlim([0 dom_size(1)])
ylim([0 dom_size(3)])

% YZ tile
nexttile
viscircles([par.r(ind_trg,2),par.r(ind_trg,3)],(par.d(ind_trg))./2, ...
    'EnhanceVisibility',false,'Color',cm1,'LineWidth',0.5);
hold on
viscircles([par.r(ind_trg,2),par.r(ind_trg,3)],...
    (coef_trg.*(par.d(ind_trg))./2)+(max(par.d)/2),...
    'EnhanceVisibility',false,'Color','k','LineStyle','--','LineWidth',0.5);
viscircles([par.r(ind_nn,2),par.r(ind_nn,3)],(par.d(ind_nn))./2, ...
    'EnhanceVisibility',false,'Color',cm2,'LineWidth',0.5);
viscircles([par.r(ind_rest,2),par.r(ind_rest,3)],(par.d(ind_rest))./2, ...
    'EnhanceVisibility',false,'Color',cm3,'LineWidth',0.5);
hold off
axis equal
title('yz view')
xlabel('y (m)')
ylabel('z (m)')
xlim([0 dom_size(2)])
ylim([0 dom_size(3)])

% XYZ tile
nexttile
p1_xyz = scatter3(par.r(ind_trg,1),par.r(ind_trg,2),par.r(ind_trg,3),...
    ((par.d(ind_trg))./2).*2e8,cm1,'filled');
hold on
p2_xyz = scatter3(par.r(ind_trg,1),par.r(ind_trg,2),par.r(ind_trg,3),...
    ((coef_trg.*(par.d(ind_trg))./2)+(max(par.d)/2)).*2e8,...
    'k','filled','MarkerFaceAlpha',.4);
p3_xyz = scatter3(par.r(ind_nn,1),par.r(ind_nn,2),par.r(ind_nn,3),...
    ((par.d(ind_nn))./2).*2e8,cm2,'filled');
p4_xyz = scatter3(par.r(ind_rest,1),par.r(ind_rest,2),par.r(ind_rest,3),...
    ((par.d(ind_rest))./2).*2e8,cm3,'filled');
hold off
axis equal
title('xyz view')
xlabel('x (m)')
ylabel('y (m)')
zlabel('z (m)')
lgd2 = legend([p1_xyz p2_xyz p3_xyz p4_xyz],{'Target','Limit','Neighbor','Rest'});
lgd2.Layout.Tile = 'east';
xlim([0 dom_size(1)])
ylim([0 dom_size(2)])
zlim([0 dom_size(3)])

end

