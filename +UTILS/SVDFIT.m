function [outputArg1,outputArg2] = SVDFIT(inputArg1,inputArg2)
%SVDFIT Summary of this function goes here
%   Detailed explanation goes here
outputArg1 = inputArg1;
outputArg2 = inputArg2;

% Assuming `log_x` and `log_y` are log-transformed data points

% Step 1: Perform SVD on the data matrix
A = [log_x, log_y];
[~, S, V] = svd(A, 'econ');
direction_vector = V(:, end);  % Direction vector with smallest singular value
a_tls = -direction_vector(1) / direction_vector(2);  % TLS slope
b_tls = mean(log_y) - a_tls * mean(log_x);           % TLS intercept

% Step 2: Calculate orthogonal residuals and variance
residuals = (log_y - (a_tls * log_x + b_tls)) ./ sqrt(1 + a_tls^2);
residual_variance = var(residuals);

% Step 3: Estimate standard errors for slope and intercept
n = length(log_x); % Number of data points
x_mean = mean(log_x);
Sxx = sum((log_x - x_mean).^2);
slope_se = sqrt(residual_variance / Sxx);  % Standard error of the slope
intercept_se = sqrt(residual_variance * (1/n + x_mean^2 / Sxx));  % Intercept SE

% Step 4: Calculate the 95% confidence intervals
t_critical = tinv(0.975, n - 2); % Critical t-value for 95% CI with (n-2) degrees of freedom

slope_ci = [a_tls - t_critical * slope_se, a_tls + t_critical * slope_se];
intercept_ci = [b_tls - t_critical * intercept_se, b_tls + t_critical * intercept_se];

% Display results
fprintf('TLS Slope: %.4f with 95%% CI: [%.4f, %.4f]\n', a_tls, slope_ci(1), slope_ci(2));
fprintf('TLS Intercept: %.4f with 95%% CI: [%.4f, %.4f]\n', b_tls, intercept_ci(1), intercept_ci(2));

% Optional: Plot data with confidence intervals
scatter(log_x, log_y, '.');
hold on;
plot(log_x, a_tls * log_x + b_tls, 'r', 'LineWidth', 1.5);
plot(log_x, (a_tls - t_critical * slope_se) * log_x + (b_tls - t_critical * intercept_se), 'r--');
plot(log_x, (a_tls + t_critical * slope_se) * log_x + (b_tls + t_critical * intercept_se), 'r--');
xlabel('log(X)');
ylabel('log(Y)');
title('TLS Fit with 95% Confidence Interval');
legend('Data', 'TLS Fit', '95% Confidence Interval');
hold off;

end

