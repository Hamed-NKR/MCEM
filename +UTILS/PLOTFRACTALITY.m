function [df1, kf1, h_fract] = PLOTFRACTALITY(parsdata, kk)
% "PLOTFRACTALITY" can plot the instantaneous as well as compiled...
%   ...fractal properties of an aggregate population.
% ----------------------------------------------------------------------- %
%
% Input:
%   parsdata: Stored particle information over time
%   kk: Marching indices of stored data to be plotted
% ----------------------------------------------------------------------- %
%
% Outputs:
%   df_compiled: Fractal dimension
%   kf_compiled: Fractal prefactor
%   h_fract: The output figure handle
% ----------------------------------------------------------------------- %

% Setting default indices of data to be plotted
if ~exist('kk', 'var'); kk = []; end
if isempty(kk)
    kk = unique(round(1 : (length(parsdata.ii) - 1) / 10 :...
        length(parsdata.ii)));
else
    kk = unique(kk);
end

% Initializing the results plot
figure;
h_fract = gcf;
if ~all(h_fract.Position == [0, 0, 2000, 892.1])
    h_fract.Position = [0, 0, 2000, 892.1];
end
set(h_fract, 'color', 'white');

% Changes of number and collision frequency over time
subplot(1,2,1)
yyaxis left
plot(parsdata.t, parsdata.df);
xlabel('t (s)')
ylabel('df (-)')
ylim([1 3])
yyaxis right
plot(parsdata.t, parsdata.kf);
ylabel('kf (-)')
ylim([1 inf])
title ('Fractal properties over time')
set(gca, 'XScale', 'log')

% Setting a color distribution for the second plot lines
cl = jet;
ii = unique(round(10 + (length(cl) - 10) .*...
    (0 : 1 / (length(kk) - 1) : 1)'));
cl = cl(ii,:);

% Plotting the fractal scattered data
subplot(1,2,2)
ms = 25; % Marker size
for i = 1 : length(kk)
    scatter(parsdata.dg_dpp{kk(i)}, parsdata.npp{kk(i)}, ms, cl(i,:),...
        'filled');
    hold on
end

legtxt = "t = " + num2str(parsdata.t(kk), '%1.1e') + " (s)"; % Time legend

% Compiling number distribution and size ratio across aggregates
dg_dpp = cat(1, parsdata.dg_dpp{kk});
npp = cat(1, parsdata.npp{kk});

%%% Data fits
% set(gca, 'XScale', 'log')
% set(gca, 'YScale', 'log') % Logarithmic axes scales

% #0. DLCA
df0 = 1.78;
kf0 = 1.30; % DLCA properties from Sorensen (2011)
xfit = min(dg_dpp) : range(dg_dpp) / (10 * numel(dg_dpp)) : max(dg_dpp);
yfit0 = kf0 .* xfit .^ df0;

anotxt0 = "df0 = " + num2str(df0, '%.2f') +...
    ", kf0 = " +num2str(kf0, '%.2f');
cl0 = [0.4660 0.6740 0.1880];
plot(xfit, yfit0, 'Color', cl0, 'LineStyle', '-', 'LineWidth', 2);
    % Plotting the fit
annotation('textbox', [0.8 0.1 1 0.08], 'String', anotxt0,...
    'Color', cl0, 'EdgeColor', cl0, 'FitBoxToText', 'on');
        % Fractal properties annotated

% #1. Polynomial
fit1 = polyfit(log(dg_dpp), log(npp), 1); % Fitting log(y) = b*log(x)+log(a) 
df1 = fit1(1); % Compiled fractal dimension
kf1 = exp(fit1(2)); % Compiled fractal prefactor
yfit1 = kf1 .* xfit .^ df1;

anotxt1 = "df1 = " + num2str(df1, '%.2f') +...
    ", kf1 = " +num2str(kf1, '%.2f');
cl1 = [0.6350 0.0780 0.1840];
plot(xfit, yfit1, 'Color', cl1, 'LineStyle', '--');
annotation('textbox', [0.8 0.12 1 0.1], 'String', anotxt1,...
    'Color', cl1, 'EdgeColor', cl1, 'FitBoxToText', 'on');

% #2. Power series (without intercept)
fit2 = fit(dg_dpp, npp, 'power1'); % fitting y = a*x^b
df2 = fit2.b;
kf2 = fit2.a;
y_fit2 = kf2 .* xfit .^ df2;

anotxt2 = "df2 = " + num2str(df0, '%.2f') +...
    ", kf2 = " +num2str(kf0, '%.2f');
cl2 = [0.9290 0.6940 0.1250];
plot(xfit, y_fit2, 'Color', cl2, 'LineStyle', '-.');
annotation('textbox', [0.8 0.14 1 0.12], 'String', anotxt2,...
    'Color', cl2, 'EdgeColor', cl2, 'FitBoxToText', 'on');

% #3. Power series (with intercept)
fit3 = fit(dg_dpp, npp, 'power2'); % fitting y = a*x^b + c
df3 = fit3.b;
kf3 = fit3.a;
y_fit3 = kf3 .* xfit .^ df3;

cl3 = [0 0.4470 0.7410];
anotxt3 = "df3 = " + num2str(df0, '%.2f') +...
    ", kf3 = " +num2str(kf0, '%.2f');
plot(xfit, y_fit3, 'Color', cl3, 'LineStyle', '-.');
annotation('textbox', [0.8 0.16 1 0.14], 'String', anotxt3,...
    'Color', cl3, 'EdgeColor', cl3, 'FitBoxToText', 'on');

% #4. Linear regression
fit4 = fitlm(table(log(dg_dpp), log(npp)), 'linear'); % Fitting...
    % ...log(y) = b*log(x)+log(a)
df4 = fit4.Coefficients.Estimate(2);
kf4 = exp(fit4.Coefficients.Estimate(1));
y_fit4 = kf4 .* xfit .^ df4;

cl4 = [0.4940 0.1840 0.5560];
anotxt4 = "df4 = " + num2str(df0, '%.2f') +...
    ", kf4 = " +num2str(kf0, '%.2f');
plot(xfit, y_fit4, 'Color', cl4, 'LineStyle', ':');
annotation('textbox', [0.8 0.18 1 0.16], 'String', anotxt4,...
    'Color', cl4, 'EdgeColor', cl4, 'FitBoxToText', 'on');

legtxt = [legtxt; 'DLCA benchmark'; 'Linear polynomial';...
    'Power law 1'; 'Power law 2'; 'Linear regression'];

title('Time- & aggregate-ensembled fractal properties')
xlabel('dg / dpp,g (-)')
ylabel('npp (#)')
legend(legtxt, 'Location', 'northwest')
axis padded

disp(' ')
disp('Compiled average fractal properties:')
fprintf('df_ens = %.2f \n', df1)
fprintf('kf_ens = %.2f \n', kf1)

if nargout < 3
    clear h_fract;  % Deleting figure handle if not requested as an output
end

end

