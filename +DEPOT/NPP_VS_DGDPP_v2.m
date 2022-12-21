function h = NPP_VS_DGDPP_v2(parsdata, opts)
% "PA_VS_NPP" plots the number of primaries vs. gyration diameter over...
%   ...primary particle diameter for two sets of hybrid and non-hybrid...
%   ...aggregates and compares them with benchmark values.
% ----------------------------------------------------------------------- %
% 
% Inputs:
%   parsdata: a cell array of structures containing temporal aggregate info
%   opts: plotting options
% ----------------------------------------------------------------------- %
%
% Outputs:
%   h: output figure handle
% ----------------------------------------------------------------------- %

% initialize figure 
figure;
h = gcf;
h.Position = [0, 0, 800, 800];
set(h, 'color', 'white');

% initialize plotting options setting
if ~exist('opts', 'var') 
    opts = struct();
end

% field for data fitting and getting fractal properties
if ~isfield(opts, 'fit')
    opts.fit = [];
end
opts_fit = opts.fit;
if isempty(opts_fit)
    opts_fit = 'off'; % default not to fit the data
end

% field for dataset to be selected for fitting
if ~isfield(opts, 'dat')
    opts.dat = [];
end
opts_dat = opts.dat;
if isempty(opts_dat)
    opts_dat = 0; % default to ensemble fit(throughout all timesteps)
end

% field for the type of fitting to be performed on the ensemble data
if ~isfield(opts, 'type')
    opts.type = [];
end
opts_type = opts.type;
if isempty(opts_type)
    opts_type = 'none'; % default not to filter the data for fitting
end

% field for post-filtering number of aggs
if ~isfield(opts, 'n_agg')
    opts.n_agg = [];
end
opts_nagg = opts.n_agg;

n_dat = length(parsdata); % number of hybrid datasets to be plotted

% data extents
kk_max = 0.5; 
kk_min = 0.02;

r_t = (kk_min / kk_max)^(1 / (n_dat - 2));
t_id = kk_max * ones(n_dat - 1,1);
for i = 2 : (n_dat - 1)
    t_id(i) = t_id(i) * r_t^(i-1);
end
t_id = [1; t_id]; % time identifier of data to be plotted 

if ismember(opts_fit, {'ON', 'On', 'on'})
    p = cell(n_dat + 3,1); % initialize the plot cell array
    legtxt = cell(n_dat + 3,1); % placeholder for legends
else
    p = cell(n_dat + 1,1);
    legtxt = cell(n_dat + 1,1);
end

% set colormap
mc = colormap(hot);
ii = round(1 + (length(mc) - 1) .* (0.05 : 0.8 / (n_dat - 1) : 0.85)');
mc = mc(ii,:);
% mc = flip(mc,1);

ms = [25, 25, 25, 55, 35, 60]; % Marker sizes
mt = {'o', '^', 'v', 's', 'd', 'p'}; % Marker types

% plot the Meakin et al. (1989) recommended correlation
rd0 = log10(logspace(0, 1e3, 1e4));
n0 = 1.3 * rd0.^1.78;
p{end} = plot(rd0, n0, 'Color', [0.4940 0.1840 0.5560], 'LineStyle', '-.',...
    'LineWidth', 2.5);
legtxt{end} = strcat(' \itd_{f_0}\rm = 1.78 & \itk_{f_0}\rm = 1.30 (Sorensen, 2011)');
% legtxt{end} = strcat(' \itd_{f_0}\rm = 1.78 & \itk_{f_0}\rm = 1.30',...
%     string(newline), '(Sorensen, 2011)');
hold on

pltvar = false; % decision variable for plotting

% plot temporal aggregate area vs. number data
for i = 1 : n_dat
    p{i} = scatter(parsdata(i).dg ./ parsdata(i).dpp_g(:,1),...
        parsdata(i).npp, ms(i), mc(i,:), mt{i}, 'LineWidth', 0.1);
    
    legtxt{i} = strcat(' \itn_{agg}\rm/\itn_{agg_0}\rm =',...
        {' '}, num2str(t_id(i), '%.2f'));
    
    if ismember(opts_fit, {'ON', 'On', 'on'})
        if i == opts_dat % assign the dataset that needs to be fitted
            fit = fitlm(table(log(parsdata(i).dg ./ parsdata(i).dpp_g(:,1)),...
                log(parsdata(i).npp)), 'linear'); % Fit a linear regression...
                    % ...model to the population
            
            df = fit.Coefficients.Estimate(2); % fractal dimension
            kf = exp(fit.Coefficients.Estimate(1)); % prefactor            
            
            legtxt{end - 2} = strcat(' \itd_{f_{', num2str(t_id(i), '%.2f'),...
                '}}\rm =', {' '}, num2str(df, '%.2f'), ' & \itk_{f_{',...
                num2str(t_id(i), '%.2f'), '}}\rm =', {' '}, num2str(kf, '%.2f'));
            
            pltvar = ~pltvar; % update plotting status
            
        elseif ~pltvar && (opts_dat <= 0)
            
            % concatinate the data
            dpp_ens = cat(1, parsdata.dpp_g);
            dpp_ens = dpp_ens(:,1);
            dgg_dpp_ens = cat(1, parsdata.dg) ./ dpp_ens;
            npp_ens = cat(1, parsdata.npp);
            da_ens = cat(1, parsdata.da);
                        
            % uniform randam sampling of ensemble data
            if ismember(opts_type, {'UNIFORM', 'Uniform', 'uniform'})
                if isempty(opts_nagg); opts_nagg = length(npp_ens); end 
                ind = randperm(length(npp_ens), opts_nagg); % index of particles selected
                
                % random aggregate sampling
                dgg_dpp_fit = dgg_dpp_ens(ind);
                npp_ens_fit = npp_ens(ind);
                
            % lognormal sampling
            elseif ismember(opts_type, {'LOGNORMAL', 'Lognormal', 'lognormal'})
                if isempty(opts_nagg); opts2.nfit = 'on'; end 
                [~, ind] = TRANSP.LNSAMPLING(da_ens, median(da_ens), 1.4,...
                    opts_nagg, 2^4, [], opts2);
                
                ind = ind(~cellfun('isempty', ind));
                ind = cat(1, ind{:});
                
                % filter aggregates
                dgg_dpp_fit = dgg_dpp_ens(ind);
                npp_ens_fit = npp_ens(ind);
                
            % pass a fit from the middle of data
            elseif ismember(opts_type, {'MEDIAN', 'Median', 'median'})
                n_bin = 2^7; % default number of bins
                
                % set the bin locations
                dd = [min(dgg_dpp_ens), max(dgg_dpp_ens)];
                rr = (dd(2) / dd(1))^(1 / n_bin);
                d_bin = dd(1) * ones(n_bin + 1, 1);
                dd_c = zeros(n_bin,1); % bin centers
                npp_c = zeros(n_bin,1); % fit array
                
                for ii = 1 : n_bin
                    d_bin(ii+1) = d_bin(ii+1) * rr^(ii); % bin upper side 
                    dd_c(ii) = sqrt(d_bin(ii) * d_bin(ii+1)); % bin center
                    iii = (dgg_dpp_ens >= d_bin(ii)) &...
                        (dgg_dpp_ens < d_bin(ii+1)); % particles that fall in the bin
                    npp_c(ii) = median(npp_ens(iii)); % center of dependent variable in current bin
                end
                
                % assign the medians
                dgg_dpp_fit = dd_c;
                npp_ens_fit = npp_c;
                
            else
                % fit type selected as 'none':
                dgg_dpp_fit = dgg_dpp_ens;
                npp_ens_fit = npp_ens;
            end
            
            if ~ismember(opts_type, {'WEIGHTED', 'Weighted', 'weighted'})
                
                fit = fitlm(table(log(dgg_dpp_fit), log(npp_ens_fit)), 'linear'); % fit a linear regression

                df = fit.Coefficients.Estimate(2); % ensemble fractal dimension
                kf = exp(fit.Coefficients.Estimate(1)); % ~ prefactor

                legtxt{end - 2} = strcat(' \itd_{f_{ens}}\rm =', {' '}, num2str(df, '%.2f'),...
                    ' & \itk_{f_{ens}}\rm =', {' '}, num2str(kf, '%.2f')); % legend for ensemble fit
                
                ci = coefCI(fit); % 95% confidence intervals for fractal properties
                ci_df = ci(2,:);
                ci_kf = exp(ci(1,:));

                legtxt{end - 1} = strcat(' CI_{\itd_f}\rm = (', num2str(ci_df(1), '%.2f'),...
                    {','}, num2str(ci_df(2), '%.2f'), ') & CI_{\itk_f}\rm = (',...
                    num2str(ci_kf(1), '%.2f'), {','}, num2str(ci_kf(2), '%.2f'), {')'}); % legend for CI of fit            

            else
                [fr_w, std_w] = lscov([ones(length(dgg_dpp_ens), 1), log(dgg_dpp_ens)],...
                    log(npp_ens), sqrt(npp_ens)); % weighted fitting outputs
                
                % extract fractal properties
                df = fr_w(2);
                kf = exp(fr_w(1));
                
                legtxt{end - 2} = strcat(' \itd_{f_{ens}}\rm =', {' '}, num2str(df, '%.2f'),...
                    ' & \itk_{f_{ens}}\rm =', {' '}, num2str(kf, '%.2f'));                
                
                % convert std to CI
                dci_df = 1.96 * std_w(2) / sqrt(length(dgg_dpp_ens));
                dci_kf = exp(1.96 * std_w(1) / sqrt(length(dgg_dpp_ens)));
                
                legtxt{end - 1} = strcat(' CI_{\itd_f}\rm = ±', {' '}, num2str(dci_df, '%.2e'),...
                    ' & CI_{g,\itk_f}\rm = ±', {' '}, num2str(dci_kf, '%.4f')); % CI legend for weighted fit
                
                ci_df = [df - dci_df, df + dci_df];
                ci_kf = [kf / dci_kf, kf * dci_kf];
                
            end
            
            pltvar = ~pltvar; % skip plotting for next iterations
        end
                
        if pltvar
            % generate the fit data
            npp = kf * rd0.^df;
            ci_npp = [ci_kf(1) * rd0.^ci_df(1); ci_kf(2) * rd0.^ci_df(2)];
            
            % plot the main fit and CI bounds
            p{end - 2} = plot(rd0, npp, 'Color', [0.6, 0.6, 0.6],...
                'LineStyle', '-', 'LineWidth', 2);
            p{end - 1} = plot(rd0, ci_npp(1,:), rd0, ci_npp(2,:),...
                'Color', [0.6, 0.6, 0.6], 'LineStyle', ':', 'LineWidth', 1.5);
        end
    end
end

% set axes
box on
set(gca, 'FontName', 'Calibri Light', 'FontSize', 12, 'TickLength', [0.005 0.005])
xlabel('{\itd_g} / {\itd_{pp,g}} [-]', 'FontName', 'Calibri', 'FontSize', 16)
set(gca, 'XScale', 'log')
xlim([2,3e2])
ylabel('{\itn_{pp}} [-]', 'FontName', 'Calibri', 'FontSize', 16)
ylim([1,2e4])
set(gca, 'YScale', 'log')
legend(cat(2, [p{1 : end - 2}, p{end - 1}(1), p{end}])', cat(2, legtxt{:}),...
    'FontName', 'Calibri Light', 'FontSize', 12, 'Location', 'southeast')

end

