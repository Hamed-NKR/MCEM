function [df_compiled, kf_compiled, h_fract] = PLOTFRACTALITY(parsdata, kk)
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
cl = autumn;
ii = unique(round(5 + (length(cl) - 5) .*...
    (0 : 1 / (length(kk) - 1) : 1)'));
cl = cl(ii,:);
ntot = cat(1, parsdata.ntot(kk));
cl = repelem(cl, ntot, ones(3,1));

% Compiling number distribution and size ratio across aggregates
dg_dpp = cat(1, parsdata.dg_dpp{kk});
npp = cat(1, parsdata.npp{kk});

% Plotting the fractal scattered data
subplot(1,2,2)
ms = 25; % Marker size
scatter(dg_dpp, npp, ms, cl, 'filled');
hold on

pfit = polyfit(log(dg_dpp), log(npp), 1); % Fitting a power function...
    % ...to the number vs. size ratio variations
df_compiled = pfit(1); % Compiled fractal dimension
kf_compiled = exp(pfit(2)); % Compiled fractal prefactor

% Plotting the curve fit
x_fit = min(dg_dpp) : range(dg_dpp) / (10 * numel(dg_dpp)) : max(dg_dpp);
y_fit = kf_compiled .* x_fit .^ df_compiled;
plot(x_fit, y_fit, 'Color', 'r', 'LineStyle', '--');
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
title('Time- & aggregate-ensembled fractal properties')
xlabel('dg / dpp,g (-)')
ylabel('npp (#)')
axis padded
txt_compiled = "df(ens) = " + num2str(df_compiled, '%.1f') + ", "...
    + newline + " kf(ens) = " +num2str(df_compiled, '%.1f') + " ";
        % Overal results to be annotated on the plot
annotation('textarrow', [0.65 0.75], [0.8 0.6], 'String', txt_compiled);

disp(' ')
disp('Compiled average fractal properties:')
fprintf('df_ens = %.2f \n', df_compiled)
fprintf('kf_ens = %.2f \n', kf_compiled)

if nargout < 3
    clear h_fract;  % Deleting figure handle if not requested as an output
end

end

