function h = PCF_MULTI_v2(g_multi, r_multi, pp_multi)
% PCF_MULTI gets multiple pair-correlation function (PCF) distributions...
%   ...from different scenarios and organizes them for comparison...
%   ...in a table of PCF vs. radial location plots. The plots...
%   ...are separated based on the number of primaries, polydispersity...
%   ...and hybridity of aggregates/agglomerates.
% ----------------------------------------------------------------------- %
% 
% Inputs:
%   g_multi: An m*n*l cell array containing PCF data. Each cell contains...
%       ...data from a certain trial of a certain scenario of hybridity;...
%       ...i.e, (m) denotes hybridity level, (n) denotes the controlled...
%       ...polydispersity and size, and (l) denotes whether we look at...
%       ...data from a real hybrid agglomerate or one of the repeats of...
%       ...a uniform aggregate with similarly initialized parameters...
%       ...(number of primaries, geometric standard deviation of primary
%       ...particle size).
%   r_multi: An m*n*l cell array including the raidal logarithmic points...
%       ...of calculation corresponding to the data in g_multi.
%   pp_multi: An m*n*l cell array of primary particle info matrices...
%       ...corresponding to the g_multi+r_multi data.
% ----------------------------------------------------------------------- %
% 
% Outputs:
%   h: PCF plot figure handle
% ----------------------------------------------------------------------- %

% number of rows and columns in subplots and number of lines within them
[m,n,l] = size(g_multi);

% initialize figure 
figure;
h = gcf;
h.Position = [0, 0, n * 500, m * 500];
set(h, 'color', 'white');

% initialize layout
tt = tiledlayout(m,n);
tt.TileSpacing = 'compact';
tt.Padding = 'compact';

lin_c1 = [203,24,29;
    33,113,181;
    35,139,69;
    106,81,163;
    217,71,1;
    82,82,82] / 255; % line color repository - stream 1

lin_c2 = [252,174,145;
    189,215,231;
    186,228,179;
    203,201,226;
    253,190,133;
    204,204,204] / 255; % line color repository - stream 2

lin_w = [3, 1.5]; % line width repo.
lin_t = {'-', ':'}; % line style repo.

for i = 1 : m
    for j = 1 : n
        nexttile
        
        % plot the hybrid agglomerates
        plt1 = plot(r_multi{i,j,1}, g_multi{i,j,1}, 'Color', lin_c1(i,:),...
            'LineStyle', lin_t{1}, 'LineWidth', lin_w(1));
        
        hold on
        
        % plot the equivalent uniform aggregates
        for k = 2 : l
            plt2 = plot(r_multi{i,j,k}, g_multi{i,j,k}, 'Color', lin_c2(i,:),...
                'LineStyle', lin_t{2}, 'LineWidth', lin_w(2));
        end
        
        if (i == 1) && (j == n)
            legend([plt1, plt2], {'Hybrid', 'Uniform'},...
                'interpreter', 'latex', 'FontSize', 16,...
                'Location', 'northeast');
        end
        
        % make the subtitles with properties of agglomerates
        nhyb = length(unique(pp_multi{i,j,1}(:,6)));
        sigma = UTILS.GEOSTD(pp_multi{i,j,1}(:,2));
        npp = size(pp_multi{i,j,1},1);
        ttl = strcat({'$n_\mathrm{hyb}$ ='}, {' '}, num2str(nhyb, '%d'),...
            {','}, ' $\sigma_\mathrm{agg}$ =', {' '}, num2str(sigma, '%.3f'),...
            {','}, {' $n_\mathrm{pp}$ ='}, {' '}, num2str(npp, '%d'));
        title(ttl, 'FontSize', 18, 'interpreter','latex')
        
        box on
        set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 18,...
            'TickLength', [0.02 0.02], 'XScale', 'log', 'YScale', 'log')
        xlim([0.1 120])
        ylim([5e-5 1])
        ax = gca;
        set(ax, 'xtick',[0.1, 1, 10, 100])
        set(ax, 'ytick',[1e-4, 1e-3, 1e-2, 1e-1, 1])
        if ismember(j, 2:n)
            set(ax, 'yticklabel',[])
        end
        if ismember(i, 1:(m-1))
            set(ax, 'xticklabel',[])
        end
    end
end

xlabel(tt, '$\hat{r}$ (-)', 'interpreter', 'latex', 'FontSize', 20)
ylabel(tt, '$\overline{g}$($\hat{r}$) (-)', 'interpreter', 'latex', 'FontSize', 20)

end

