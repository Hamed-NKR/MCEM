function h_nn = PLOTNN(pars, dom_size, ind_trg, coef_trg)
% "PLOTNN" highlights the nearest neighbors of a desired particles.
% ----------------------------------------------------------------------- %
% 
% Inputs:
%     dom_size: The computational domain size
%     par: The particle information structure/class
%     ind_trg: The index of target particle (the one neighbors of which...
%         ...need to be identified)
%     coef_trg: The enlargement coefficient for the size of a spherical...
%         ...barrier used to identify the neighbors
% ----------------------------------------------------------------------- %
%
% Outputs
%     h_main: Figure handle for the nearest neighbors
% ----------------------------------------------------------------------- %        

% Initializing the figure handle
hold off
h_nn = gcf;

% Clearing previous data (for animations)
all_axes_in_figure = findall(h_nn, 'type', 'axes');
n_ax = numel(all_axes_in_figure);
for i = 1 : n_ax
    cla(all_axes_in_figure(i))
end

figure(h_nn);

% Setting figure position and background
h_nn.Position = [0, 0, 2000, 892.1];
set(h_nn, 'color', 'white');

% Initialing the colors of different plot elements
cm1 = [0.2,1,0.3];  % Color of the target particle
cm2 = [1,0.3,0.2];  % color of the nearest neighbors
cm3 = [0.3,0.2,1];  % Color of the non-neighbor particles

n_trg = numel(ind_trg); % Number of target particles for neighbor checking
tiledlayout(n_trg,4); % Initializing the figure layout
[ind_nn, ind_rest] = COL.NNS(pars, ind_trg, coef_trg); % Getting the...
    % ...neighbor and non-neighbor indices

% Maximum size of the particles with respect to their center of mass
if isa(pars, 'AGG')
    dmax = zeros(length(pars), 1);
    for i = 1 : length(pars)
        dmax(i) = pars.TERRITORY(pars(i).pp);
    end
    
else
    dmax = PAR.TERRITORY(pars.pp, pars.n);
end

% Concatinating the aggregates global info
r = cat(1, pars.r);

% Looping over the target particles
for i = 1 : n_trg
    
    % Compiling primary particles across multiple aggregates
    if isa(pars, 'AGG')
        pp_trg = AGG.COMPILEPP(pars(ind_trg(i)));
        pp_nn = AGG.COMPILEPP(pars(ind_nn{i}));
        pp_rest = AGG.COMPILEPP(pars(ind_rest{i}));
        
    else
        pp_trg = cell2mat(pars.pp(ind_trg(i))); % Target particles pp info
        pp_nn = cell2mat(pars.pp(ind_nn{i})); % Nearest neighbor's pp info
        pp_rest = cell2mat(pars.pp(ind_rest{i})); % Non-neighbor's pp info
        
    end
    
    % Avoiding errors while plotting due to array emptiness
    if isempty(pp_nn)
        pp_nn = NaN(1,6);
    end
    if isempty(pp_rest)
        pp_rest = NaN(1,6);
    end
    
    % XY tile
    nexttile
    % Plotting the target particle
    p1_xy = viscircles([pp_trg(:,3), pp_trg(:,4)], pp_trg(:,2) ./ 2,...
        'EnhanceVisibility', false, 'Color', cm1, 'LineWidth', 0.5);
    hold on
    % Neigboring limit
    p2_xy = viscircles([r(ind_trg(i),1), r(ind_trg(i),2)],...
        coef_trg(i) .* dmax(ind_trg(i)) ./ 2,...
        'EnhanceVisibility', false, 'Color', 'k','LineStyle','--',...
        'LineWidth',0.5);
    % Nearest neighbors
    p3_xy = viscircles([pp_nn(:,3), pp_nn(:,4)], pp_nn(:,2) ./ 2,...
        'EnhanceVisibility',false,'Color',cm2,'LineWidth',0.5);
    % Non-neighbors
    p4_xy = viscircles([pp_rest(:,3), pp_rest(:,4)], pp_rest(:,2) ./ 2,...
        'EnhanceVisibility',false,'Color',cm3,'LineWidth',0.5);
    % Setting different plot graphics
    hold off
    axis equal
    title('xy view')
    xlabel('x (m)')
    ylabel('y (m)')
    if i == 1
        lgd1 = legend([p1_xy, p2_xy, p3_xy, p4_xy],{'Targets', 'Limits',...
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
    viscircles([r(ind_trg(i),1), r(ind_trg(i),3)],...
        coef_trg(i) .* dmax(ind_trg(i)) ./ 2,...
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
    viscircles([r(ind_trg(i),2), r(ind_trg(i),3)],...
        coef_trg(i) .* dmax(ind_trg(i)) ./ 2,...
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
    p1_xyz = UTILS.PLOTPP(pp_trg(:,3), pp_trg(:,4), pp_trg(:,5),...
        pp_trg(:,2), cm1);
    hold on
    p2_xyz = UTILS.PLOTPP(r(ind_trg(i),1), r(ind_trg(i),2),...
        r(ind_trg(i),3), coef_trg(i) .* dmax(ind_trg(i)), 'k', 0.3);
    p3_xyz = UTILS.PLOTPP(pp_nn(:,3), pp_nn(:,4), pp_nn(:,5),...
        pp_nn(:,2), cm2);
    p4_xyz = UTILS.PLOTPP(pp_rest(:,3), pp_rest(:,4), pp_rest(:,5),...
        pp_rest(:,2), cm3);
    hold off
    axis equal
    title('xyz view')
    xlabel('x (m)')
    ylabel('y (m)')
    zlabel('z (m)')
    if i == 1
        lgd2 = legend([p1_xyz, p2_xyz, p3_xyz, p4_xyz],{'Targets',...
            'Limits', 'Neighbors', 'Rest'});
        lgd2.Layout.Tile = 'south';
    end
    xlim([0 dom_size(1)])
    ylim([0 dom_size(2)])
    zlim([0 dom_size(3)])
end

if nargout == 0
    clear h_nn;  % Deleting figure handle if not requested as an output
end


end

