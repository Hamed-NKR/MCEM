function [df1, kf1, h_fract] = PLOTFRACTALITY_REVISED(parsdata, kk)
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

% Instantaneous (Power series, Without intercept)
fit1 = fit(dg_dpp, npp, 'power1'); % fitting y = a*x^b
df1 = fit1.b;
kf1 = fit1.a;
xfit = min(dg_dpp) : range(dg_dpp) / (10 * numel(dg_dpp)) : max(dg_dpp);
y_fit1 = kf1 .* xfit .^ df1;

anotxt1 = "df = " + num2str(df1, '%.2f') +...
    ", kf = " +num2str(kf1, '%.2f');
cl1 = [0.4660 0.6740 0.1880];
plot(xfit, y_fit1, 'Color', cl1, 'LineStyle', '--', 'LineWidth', 1.5);
annotation('textbox', [0.8 0.14 1 0.12], 'String', anotxt1,...
    'Color', cl1, 'EdgeColor', cl1, 'FitBoxToText', 'on');

% DLCA benchmark
df0 = 1.78;
kf0 = 1.30; % DLCA properties from Sorensen (2011)
yfit0 = kf0 .* xfit .^ df0;

anotxt0 = "df_0 = " + num2str(df0, '%.2f') +...
    ", kf_0 = " +num2str(kf0, '%.2f');
cl0 = [0.5 0.5 0.5];
plot(xfit, yfit0, 'Color', cl0, 'LineStyle', '-.', 'LineWidth', 2.5);
    % Plotting the fit
annotation('textbox', [0.8 0.12 1 0.1], 'String', anotxt0,...
    'Color', cl0, 'EdgeColor', cl0, 'FitBoxToText', 'on');
        % Fractal properties annotated

legtxt = [legtxt; 'Power series'; 'Sorensen (2011)'];

title('Time- & aggregate-ensembled fractal properties')
xlabel('d_g / d_p_p (-)')
ylabel('n_p_p (-)')
legend(legtxt, 'Location', 'northwest')
axis padded

disp(' ')
disp('Average fractal properties:')
fprintf('df_ens = %.2f \n', df1)
fprintf('kf_ens = %.2f \n', kf1)

if nargout < 3
    clear h_fract;  % Deleting figure handle if not requested as an output
end

end

