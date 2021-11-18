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
    i_start = min(max(10, ceil(length(parsdata.ii) / 10)),...
        parsdata.ii(end)); % Start index for data to be plotted
    kk = unique(round(i_start : (length(parsdata.ii) - i_start) / 10 :...
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
xlabel('t (s)', 'FontName', 'Times New Roman', 'FontSize', 18,...
    'FontWeight', 'bold')
ylabel('d_f (-)', 'FontName', 'Times New Roman', 'FontSize', 18,...
    'FontWeight', 'bold')
ylim([1 3])
yyaxis right
plot(parsdata.t, parsdata.kf);
ylabel('k_f (-)', 'FontName', 'Times New Roman', 'FontSize', 18,...
    'FontWeight', 'bold')
ylim([1 inf])
title ('Time variations of fractal properties',...
    'FontName', 'Times New Roman', 'FontSize', 20, 'FontWeight', 'bold')
set(gca, 'FontName', 'Times New Roman', 'FontSize', 18, 'XScale', 'log')

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

dg_dpp(npp < 5) = []; 
npp(npp < 5) = []; 

%%% Data fits
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log') % Logarithmic axes scales

% Instantaneous (Logarithmic, Linear regression, Without intercept)
% fit1 = fit(dg_dpp, npp, 'power1'); % fitting y = a*x^b
fit1 = fitlm(table(log(dg_dpp), log(npp)), 'linear');
% df1 = fit1.b;
% kf1 = fit1.a;
df1 = [fit1.Coefficients.Estimate(2), fit1.Coefficients.SE(2)];
kf1 = [exp(fit1.Coefficients.Estimate(1)),...
    fit1.Coefficients.SE(1) * exp(fit1.Coefficients.Estimate(1))];
xfit = (min(dg_dpp) : range(dg_dpp) / (10 * numel(dg_dpp)) : max(dg_dpp))';
y_fit1 = [kf1(1) .* xfit .^ df1(1),... % the main fit
    (kf1(1) + kf1(2)) .* xfit .^ (df1(1) + df1(2)),... % upper bound
    (kf1(1) - kf1(2)) .* xfit .^ (df1(1) - df1(2))]; % lower bound

anotxt1 = "d_f = " + num2str(df1(1), '%.2f') + " +/- " +...
    num2str(df1(2), '%.2f') + newline +...
    ", k_f = " + num2str(kf1(1), '%.2f') + " +/- " +...
    num2str(kf1(2), '%.2f');
cl1 = [0.4660 0.6740 0.1880];
plot(xfit, y_fit1(:,1), 'Color', cl1, 'LineWidth', 2.5);
plot(xfit, y_fit1(:,2), xfit, y_fit1(:,3), 'LineStyle', '--',...
    'Color', cl1, 'LineWidth', 1.5);
annotation('textbox', [0.75 0.2 1 0.1], 'String', anotxt1,...
    'Color', cl1, 'EdgeColor', cl1, 'FitBoxToText', 'on',...
    'FontName', 'Times New Roman', 'FontSize', 14, 'FontWeight', 'bold');

% DLCA benchmark
df0 = 1.78;
kf0 = 1.30; % DLCA properties from Sorensen (2011)
yfit0 = kf0 .* xfit .^ df0;

anotxt0 = "d_f_,_0 = " + num2str(df0, '%.2f') +...
    ", k_f_,_0 = " +num2str(kf0, '%.2f');
cl0 = [0.5 0.5 0.5];
plot(xfit, yfit0, 'Color', cl0, 'LineStyle', '-.', 'LineWidth', 1.5);
    % Plotting the fit
annotation('textbox', [0.75 0.1 1 0.1], 'String', anotxt0,...
    'Color', cl0, 'EdgeColor', cl0, 'FitBoxToText', 'on',...
    'FontName', 'Times New Roman', 'FontSize', 14); % Fractal properties...
        % ...annotated

legtxt = [legtxt; 'Linear regression'; ''; ''; 'Sorensen (2011)'];

title('Time- & population-averaged fractal properties',...
    'FontName', 'Times New Roman', 'FontSize', 20, 'FontWeight', 'bold');
xlabel('d_g / d_p_p (-)', 'FontName', 'Times New Roman', 'FontSize', 18,...
    'FontWeight', 'bold')
ylabel('n_p_p (-)', 'FontName', 'Times New Roman', 'FontSize', 18,...
    'FontWeight', 'bold')
set(gca, 'FontName', 'Times New Roman', 'FontSize', 18)
legend(legtxt, 'Location', 'northwest', 'FontName', 'Times New Roman',...
    'FontSize', 14)
axis padded

disp(' ')
disp('Average fractal properties:')
fprintf('df_ens = %.2f +/- %.2f\n', df1(1), df1(2))
fprintf('kf_ens = %.2f +/- %.2f\n', kf1(1), kf1(2))

if nargout < 3
    clear h_fract;  % Deleting figure handle if not requested as an output
end

end

