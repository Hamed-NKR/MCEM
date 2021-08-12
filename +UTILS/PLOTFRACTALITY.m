function [df_compiled, kf_compiled, h_fract] = PLOTFRACTALITY(parsdata, kk)
% "PLOTFRACTALITY" can plot the instantaneous as well as compiled...
%   ...fractal properties of an aggregate population.
% ----------------------------------------------------------------------- %
%
% Input:
%   parsdata: Stored particle information over time
%   kk: Indices of stored data to be plotted
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
end

% Initializing the results plot
figure;
h_fract = gcf;
if ~all(h_fract.Position == [0, 0, 1000, 892.1])
    h_fract.Position = [0, 0, 1000, 892.1];
end
set(h_fract, 'color', 'white');
ms = 25; % Marker size

% Compiling number distribution and size ratio across aggregates
dg_dpp = cat(1, parsdata.dg_dpp{kk});
npp = cat(1, parsdata.npp{kk});
scatter(dg_dpp, npp, ms, 'filled');
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
title('Population-based fractal properties')
xlabel('dg / dpp,g (-)')
ylabel('npp (#)')
axis padded
txt = "kf = " + num2str(kf_compiled, '%.1f') + ", df = " +...
    num2str(df_compiled, '%.1f') + "  "; % Overal results to be...
        % ...annotated on the plot
annotation('textarrow', [0.3 0.5], [0.6 0.5], 'String', txt);

disp(' ')
disp('Compiled average fractal properties:')
fprintf('df_compiled = %.2f \n', df_compiled)
fprintf('kf_compiled = %.2f \n', kf_compiled)

if nargout < 3
    clear h_fract;  % Deleting figure handle if not requested as an output
end

end

