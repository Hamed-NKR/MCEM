function PLOTNN(dom_size, par, ind_trg, coef_trg)
% "PLOTNN" highlights the nearest neighbors of a desired particles.

% Inputs:
    % dom_size: The computational domain size
    % par: The information structure of particles population
    % ind_trg: The index of target particle (the one neighbors of which...
        % ...need to be identified)
    % coef_trg: The enlargement coefficient for the size of a spherical...
        % ...barrier used to identify the neighbors

hold off;

cm1 = [0.2,1,0.3];  % Color of the target particle
cm2 = [1,0.3,0.2];  % color of the nearest neighbors
cm3 = [0.3,0.2,1];  % Color of the non-neighbor particles

n_trg = numel(ind_trg); % Number of target particles for neighbor checking
tiledlayout(n_trg,4); % Initializing the figure layout
[ind_nn, ind_rest] = COL.NNS(par, ind_trg, coef_trg); % Getting the neighbor...
    % ...and non-neighbor indices

% Looping over the target particles
for i = 1 : n_trg
    
    pp_trg = cell2mat(par.pp(ind_trg(i))); % Target particles pp info
    pp_nn = cell2mat(par.pp(ind_nn{i})); % Nearest neighbor's pp info
    % Avoiding plotting errors due to array emptiness
    if isempty(pp_nn)
        pp_nn = NaN(1,6);
    end
    pp_rest = cell2mat(par.pp(ind_rest{i})); % Non-neighbor's pp info
    if isempty(pp_rest)
        pp_rest = NaN(1,6);
    end
    
    % XY tile
    nexttile
    p1_xy = viscircles([pp_trg(:,3), pp_trg(:,4)], pp_trg(:,2) ./ 2,...
        'EnhanceVisibility', false, 'Color', cm1, 'LineWidth', 0.5);
    hold on
    p2_xy = viscircles([par.r(ind_trg(i),1),par.r(ind_trg(i),2)],...
        coef_trg(i) .* par.d(ind_trg(i)) ./ 2,...
        'EnhanceVisibility', false, 'Color', 'k','LineStyle','--',...
        'LineWidth',0.5);
    p3_xy = viscircles([pp_nn(:,3), pp_nn(:,4)], pp_nn(:,2) ./ 2,...
        'EnhanceVisibility',false,'Color',cm2,'LineWidth',0.5);
    p4_xy = viscircles([pp_rest(:,3), pp_rest(:,4)], pp_rest(:,2) ./ 2,...
        'EnhanceVisibility',false,'Color',cm3,'LineWidth',0.5);
    hold off
    axis equal
    title('xy view')
    xlabel('x (m)')
    ylabel('y (m)')
    if i == 1
        lgd1 = legend([p1_xy p2_xy p3_xy p4_xy],{'Targets', 'Limits',...
            'Neighbors', 'Rest'});
        lgd1.Layout.Tile = 'south';
    end
    xlim([0 dom_size(1)])
    ylim([0 dom_size(2)])
    
    % XZ tile
    nexttile
    viscircles([pp_trg(:,3), pp_trg(:,5)], pp_trg(:,2) ./ 2,...
        'EnhanceVisibility', false, 'Color', cm1, 'LineWidth', 0.5);
    hold on
    viscircles([par.r(ind_trg(i),1),par.r(ind_trg(i),3)],...
        coef_trg(i) .* par.d(ind_trg(i)) ./ 2,...
        'EnhanceVisibility', false, 'Color', 'k','LineStyle','--',...
        'LineWidth',0.5);
    viscircles([pp_nn(:,3), pp_nn(:,5)], pp_nn(:,2) ./ 2,...
        'EnhanceVisibility',false,'Color',cm2,'LineWidth',0.5);
    viscircles([pp_rest(:,3), pp_rest(:,5)], pp_rest(:,2) ./ 2,...
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
    viscircles([pp_trg(:,4), pp_trg(:,5)], pp_trg(:,2) ./ 2,...
        'EnhanceVisibility', false, 'Color', cm1, 'LineWidth', 0.5);
    hold on
    viscircles([par.r(ind_trg(i),2),par.r(ind_trg(i),3)],...
        coef_trg(i) .* par.d(ind_trg(i)) ./ 2,...
        'EnhanceVisibility', false, 'Color', 'k','LineStyle','--',...
        'LineWidth',0.5);
    viscircles([pp_nn(:,4), pp_nn(:,5)], pp_nn(:,2) ./ 2,...
        'EnhanceVisibility',false,'Color',cm2,'LineWidth',0.5);
    viscircles([pp_rest(:,4), pp_rest(:,5)], pp_rest(:,2) ./ 2,...
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
    p1_xyz = scatter3(pp_trg(:,3), pp_trg(:,4), pp_trg(:,5),...
        2e9 .* pp_trg(:,2) ./ 2, cm1, 'filled');
    hold on
    p2_xyz = scatter3(par.r(ind_trg(i),1), par.r(ind_trg(i),2),...
        par.r(ind_trg(i),3), 2e9 .* coef_trg(i) .* par.d(ind_trg(i))...
         ./ 2, 'k', 'filled', 'MarkerFaceAlpha', 0.5);
    p3_xyz = scatter3(pp_nn(:,3), pp_nn(:,4), pp_nn(:,5),...
        2e9 .* pp_nn(:,2) ./ 2, cm2, 'filled');
    p4_xyz = scatter3(pp_rest(:,3), pp_rest(:,4), pp_rest(:,5),...
        2e9 .* pp_rest(:,2) ./ 2, cm3, 'filled');
    hold off
    axis equal
    title('xyz view')
    xlabel('x (m)')
    ylabel('y (m)')
    zlabel('z (m)')
    if i == 1
        lgd2 = legend([p1_xyz p2_xyz p3_xyz p4_xyz],{'Targets', 'Limits',...
            'Neighbors', 'Rest'});
        lgd2.Layout.Tile = 'south';
    end
    xlim([0 dom_size(1)])
    ylim([0 dom_size(2)])
    zlim([0 dom_size(3)])

end

end

