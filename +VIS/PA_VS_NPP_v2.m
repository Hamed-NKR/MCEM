function h = PA_VS_NPP_v2(parsdata)
% "PA_VS_NPP_v2" plots the projected area over number of promaries vs. ...
%   ...number of primaries for aggregates through their hybridization...
%   ...lifetimes. The results are also marked for number of internal...
%   ...clusters as well as standard deviation of primary particle size...
%   ...within the aggregates.
% ----------------------------------------------------------------------- %
% 
% Inputs:
%   parsdata: a cell array of structures from various simulations each...
%       ...containing temporal aggregate info
%   sigma_g: geometric standard deviation of agg population pp size
% ----------------------------------------------------------------------- %
%
% Outputs:
%   h: output figure handle
% ----------------------------------------------------------------------- %

% initialize figure 
figure;
h = gcf;
h.Position = [0, 0, 1000, 700];
set(h, 'color', 'white');

p = cell(8, 1); % initialize the plot cell array
legtxt = cell(8, 1); % placeholder for legends

% % set colormap
% mc = colormap(hot);
% ii = round(1 + (length(mc) - 1) .* (0.05 : 0.7 / (n_dat - 1) : 0.75)');
% mc = mc(ii,:);
% % mc = flip(mc,1);
% mc(6,:) = [236,230,61] / 255;

ms = [5, 8, 15, 25, 25, 40]; % Marker sizes
mt = {'o', '^', 'd', 's', 'p', '*'}; % Marker types
kk = [1, 2, 3, 5, 8, 15, inf];

% reproduce the literature benchmark correlations
% n0 = logspace(0, 10, 1e4);
r_dat0 = (1e10 / 1e0)^(1 / (1e4 - 1));
n0 = 1e0 * ones(1e4,1);
for i = 2 : 1e4
    n0(i) = n0(i) * r_dat0^(i-1);
end
cor1 = 4/pi * (0.3757 * n0 + 0.4098 * n0.^0.7689);

% set guideline properties
sigma_gl = 1 : 0.1 : 1.3;
n_gl = length(sigma_gl);
cor2 = ((0.94 + 0.03 * sigma_gl.^4.8) .* (n0.^0.46)).^2; % guidlines from...
    % ...correlation by Dastanpour & Rogak (2016)
x_gl = [5, 5.5, 5.8, 12]; % x location of guideline labels
y_gl = [0.8, 0.87, 0.91, 0.915]; % y location ~
theta_gl = [52, 52, 52, 52]; % orientation of ~
t_gl = cell(n_gl,1); % placeholder for Dastanpour & Rogak's correlation labels

% plot correlations
p{7} = plot(n0, cor1 ./ n0, 'Color', [0.4940 0.1840 0.5560],...
    'LineStyle', '-.', 'LineWidth', 4);
hold on

for ii = 1 : n_gl
    if ii == 1
        p{8} = plot(n0, cor2(:,ii) ./ n0, 'Color', [0.5 0.5 0.5],...
            'LineStyle', '--', 'LineWidth', 2.5);
        t_gl{ii} = text(x_gl(ii), y_gl(ii), strcat('$\sigma_{g,pp,ens} =$', {' '},...
             num2str(sigma_gl(ii), '%.1f'), ',', string(newline), 'from Sorensen (2011)'),...
             'interpreter', 'latex', 'FontSize', 12, 'Color', [0.1 0.1 0.1]);
        set(t_gl{ii}, 'Rotation', -theta_gl(ii));
    else
        plot(n0, cor2(:,ii) ./ n0, 'Color', [0.5 0.5 0.5],...
            'LineStyle', '--', 'LineWidth', 2.5);
        t_gl{ii} = text(x_gl(ii), y_gl(ii), strcat('$\sigma_{g,pp,ens} =$', {' '},...
             num2str(sigma_gl(ii), '%.1f')),'interpreter', 'latex',...
             'FontSize', 12, 'Color', [0.1 0.1 0.1]);
        set(t_gl{ii}, 'Rotation', -theta_gl(ii));
    end
end

legtxt{7} = 'Meakin et al. (1989)';
legtxt{8} = 'Dastanpour \& Rogak (2016)';

%%% plot concatinated aggregate area vs. number data

% combine data from different simulations (different ensemble standard...
    % ...deviations of primary particle size)
n0 = length(parsdata);
dpp = [];
da = [];
npp = [];
nhyb = [];
pp = [];
for ii = 1 : n0
    dpp = [dpp; cat(1, parsdata{ii}.dpp_g)];
    da = [da; cat(1, parsdata{ii}.da)];
    npp = [npp; cat(1, parsdata{ii}.npp)];
    pp = [pp; cat(1, parsdata{ii}.pp)];
    if isempty(parsdata{ii}(1).n_hyb)
        parsdata{ii}(1).n_hyb = ones(length(parsdata{ii}(1).npp),1);
    end
    nhyb = [nhyb; cat(1, parsdata{ii}.n_hyb)];
end

% combine properties from within and between various simulations
dpp_ens = cat(1, dpp);
da_ens = cat(1, da);
npp_ens = cat(1, npp);
nhyb_ens = cat(1, nhyb);
pp_ens = cat(1, pp);
nagg_ens = length(npp_ens);

% find and remove duplicates
ij = nchoosek(1 : nagg_ens, 2);
ind_flt = zeros(nagg_ens,1);
for i = 1 : length(ij)
    if isequal(sort(unique(pp_ens{ij(i,1)}(:,1))),...
            sort(unique(pp_ens{ij(i,2)}(:,1))))
        ind_flt(ij(i,2)) = 1;
    end
end
% % ind_flt = zeros(nagg_ens,1);
% % ind_flt = logical(ind_flt);
% dpp_flt = dpp_ens(~ind_flt,1);
sigmapp_flt = dpp_ens(~ind_flt,2);
da_flt = da_ens(~ind_flt);
npp_flt = npp_ens(~ind_flt);
nhyb_flt = nhyb_ens(~ind_flt);
pp_flt = pp_ens(~ind_flt);
nagg_ens = length(npp_flt);

% get mean primary particle area for normalization
dpp0 = zeros(nagg_ens, 1);
for k = 1 : nagg_ens
    dpp0(k) = sqrt(sum(pp_flt{k}(:,2)).^2) / npp_flt(k);
end

iii = cell(6,1); % index sorting placeholder based on number of internal clusters

for i = 1 : 6
    iii{i} = (nhyb_flt >= kk(i)) & (nhyb_flt < kk(i+1));
    
    p{i} = scatter(npp_flt(iii{i}), ((da_flt(iii{i}) ./ dpp0(iii{i})).^2)...
        ./ npp_flt(iii{i}), ms(i), sigmapp_flt(iii{i}), mt{i},...
        'LineWidth', 1);

    % make legends for pannel 2
    switch i
        case {1,2}
            legtxt{i} = strcat('$n_{hyb}$ =',...
                {' '}, num2str(kk(i), '%d'));
        case {3,4,5}
            legtxt{i} = strcat(num2str(kk(i), '%d'), {' '},...
                '$\leq n_{hyb} <$', {' '}, num2str(kk(i+1), '%d'));
        otherwise
            legtxt{i} = strcat('$n_{hyb} \geq$', {' '},...
                num2str(kk(i), '%d'));
    end
end

% dummy point to adjust colorbar
scatter(1e4, 0.9, 10, 1.6,'*')

% set subplots' properties
box on
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 18,...
    'TickLength', [0.02 0.02], 'XScale', 'log', 'YScale', 'log')
xlim([4, 1e4])
ylim([0.54, 0.92])

cb = colorbar;
colormap turbo
cb.Label.String = '$\sigma_{g,pp,agg}$ [-]';
cb.Label.Interpreter  = 'latex';
% cb.Label.Rotation = 360;
cb.TickLabelInterpreter  = 'latex';
cb.FontSize = 16;
cb.Label.FontSize = 20;
cb.Limits = [1.0 1.6];
cb.Ticks = 1.0:0.1:1.6;
cb.TickLength = 0.02;
cb.Location = 'northoutside';
% cbpos = get(cb, 'Position');
% cb.Label.Position = [cbpos(1) - 0.75 , cbpos(2) + 1.705];

% dummy plots for the legend appearance purposes
pnot = cell(6,1);
nppnot = 1e-3 * (1 : 6);
Abarnot = 1e-3 * (1 : 6);
for i = 1 : 6
    pnot{i} = scatter(nppnot(i), Abarnot(i), ms(i), [0.1 0.1 0.1], mt{i});
end

% set general plot's properties
xlabel('$n_{pp}$ [-]', 'interpreter', 'latex', 'FontSize', 20)
ylabel('$\hat{\overline{A}}_{agg}$ [-]', 'interpreter', 'latex', 'FontSize', 20)
lgd = legend(cat(2, [pnot{:}, p{7:8,1}])', cat(2, legtxt{:})',...
    'interpreter', 'latex', 'FontSize', 16);
lgd.Location = 'northeastoutside';

end
