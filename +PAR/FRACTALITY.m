function [df, kf, h_fract] = FRACTALITY(pars)
% "FRACTALITY" estimates the fractal properties of an aggregate population.
% ----------------------------------------------------------------------- %
%
% Inputs:
%   pars: Particle information structure/class
% ----------------------------------------------------------------------- %
%
% Output:
%   df: Fractal dimension
%   kf: Fractal prefactor
%   h_fract: the output figure handle
% ----------------------------------------------------------------------- %

% Initializing the results plot
figure;
h_fract = gcf;
if ~all(h_fract.Position == [0, 0, 1000, 892.1])
    h_fract.Position = [0, 0, 1000, 892.1];
end
set(h_fract, 'color', 'white');
ms = 25; % Marker size

% Compiling properties across aggregates
dg = cat(1, pars.dg);
dpp = cat(1, pars.dpp);
n = cat(1, pars.n);

c_d = dg ./ dpp(:,1); % Aggregate to primary particles size ratio
scatter(c_d, n, ms, 'filled');
hold on

pfit = polyfit(log(c_d), log(n), 1); % Fitting a power function...
    % ...to the number vs. size ratio variations
df = pfit(1); % Fractal dimension
kf = exp(pfit(2)); % Fractal prefactor

% Plotting the curve fit
x_fit = min(c_d) : range(c_d) / (10 * numel(c_d)) : max(c_d);
y_fit = kf .* x_fit .^ df;
plot(x_fit, y_fit, 'Color', 'r', 'LineStyle', '--');
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
title('Population-based fractal properties')
xlabel('dg / dpp,g')
ylabel('npp')
axis padded
txt = "kf = " + num2str(kf, '%.1f') + ", df = " + num2str(df, '%.1f') +...
    "  "; % Overal results to be annotated on the plot
annotation('textarrow', [0.3 0.5], [0.6 0.5], 'String', txt);

disp(' ')
disp('Averaged fractal properties:')
fprintf('df = %.1f \n', df)
fprintf('kf = %.1f \n', kf)

if nargout < 3
    clear h_fract;  % Deleting figure handle if not requested as an output
end

end

